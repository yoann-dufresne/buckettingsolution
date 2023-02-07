#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


class FastaStream
{
public:
	ifstream fs;

	FastaStream (string filename) {
		this->fs.open(filename);
		if (not this->fs.is_open()) {
			cerr << "Impossible to open " << filename << endl;
			exit(1);
		}

		// Read the first header
		char c = '\0';
		int i=0;
		while (c != '\n') {
			this->fs.get(c);
		}
	}

	~FastaStream () {
		this->fs.close();
	}

	bool has_next() {
		return ! this->fs.eof();
	}

	char next_char() {
		while (true) {
			if (! this->has_next())
				return '\0';

			char c = '\0'; bool done = false;
			this->fs.get(c);

			switch (c) {
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
				return c;
			case '>':
				while (c != '\n' and this->has_next())
					this->fs.get(c);
				return '\0';
			}
		}
	}
};


class KmerStream
{
public:
	size_t k;
	uint64_t current_kmer;
	// Mask to remove the 2 most significant bits
	uint64_t k1_mask;

	KmerStream(size_t k) : k(k), current_kmer(0)
						 , k1_mask((1ul << (2 * (k-1))) - 1)
	{};

	uint64_t next_kmer (char nucl) {
		// Remove the leftmost nucleotide
		this->current_kmer &= this->k1_mask;
		// Shift the kmer to the left creating a "hole" on the right
		this->current_kmer <<= 2;
		// Add the new nucleotide
		switch (nucl) {
		case 'A':
		case 'a':
			this->current_kmer += 0;
			break;
		case 'C':
		case 'c':
			this->current_kmer += 1;
			break;
		case 'G':
		case 'g':
			this->current_kmer += 2;
			break;
		case 'T':
		case 't':
			this->current_kmer += 3;
			break;
		}

		return this->current_kmer;
	}
};


#include <algorithm>
string uint2kmer (uint64_t val, size_t k)
{
	stringstream ss;
	for (size_t i=0 ; i<k ; i++)
	{
		char c;
		switch (val & 0b11) {
		case 0: c = 'A'; break;
		case 1: c = 'C'; break;
		case 2: c = 'G'; break;
		case 3: c = 'T'; break;
		}
		ss << c;
		val >>= 2;
	}

	string kmer = ss.str();
	reverse(kmer.begin(), kmer.end());
	return kmer;
}


#include <sys/stat.h>
class SkmerSaver
{
public:
	string outdir;
	unordered_map<uint64_t, stringstream> file_buffers;

	SkmerSaver(string outdir) : outdir(outdir)
	{
		// Verify the existance of the out directory
		struct stat sb;

		if (stat(outdir.c_str(), &sb) != 0)
		{
			cerr << "Invalid output folder for " << outdir << endl;
			cerr << "Maybe the folder does not exist ?" << endl;
			exit(1);
		}

		if (outdir[outdir.length() - 1] != '/')
			this->outdir = outdir + "/";
	};

	~SkmerSaver()
	{
		for (auto & pair : this->file_buffers)
		{
			const uint64_t minimizer = pair.first;
			ofstream mini_file(this->outdir + to_string(minimizer) + ".txt");
			string val(pair.second.str());
			mini_file.write(val.c_str(), val.length());
			mini_file.close();
		}
	}

	/** Uses the circular buffer from the inputs to save a superkmer in the right buffer.
	 **/
	void save_skmer (const uint64_t minimizer, const char * buffer, const size_t buff_size, size_t start, size_t stop)
	{
		if (stop == start)
			return;

		// cout << "save " << start << " " << stop << endl;

		// Create a new string stream if not already present
		if (this->file_buffers.find(minimizer) == this->file_buffers.end())
			this->file_buffers[minimizer] = stringstream();

		// Construct the superkmer in the right buffer
		for (size_t i=start ; i<=stop ; i++)
			this->file_buffers[minimizer] << buffer[i % buff_size];
		this->file_buffers[minimizer] << endl;
	};
};


void compute_skmers (string filename, string outdir, const size_t k, const size_t m)
{
	FastaStream fasta(filename);
	KmerStream stream(m);
	SkmerSaver * saver = new SkmerSaver(outdir);

	// Number of nucleotides to read before the first kmer
	size_t position = 0;
	const size_t nb_mmers_in_kmer = k - m + 1;
	const size_t max_skmer_size = 2 * k - m;

	// Loop variables
	uint64_t minimizer = 1ul << (2 * m);
	size_t minimizer_pos = 0;
	uint64_t * candidates = new uint64_t[nb_mmers_in_kmer];
	char * skmer_buffer = new char[max_skmer_size];
	size_t skmer_start = 0;

	char c = fasta.next_char();
	while (fasta.has_next())
	{
		// New sequence
		if (c == '\0')
		{
			position = 0;
			minimizer = 1ul << (2 * m);
			skmer_start = 0;
		}
		// New char in current sequence
		else
		{
			uint64_t candidate = stream.next_kmer(c);
			
			// cout << "char " << c << " [";
			// for (size_t i=0 ; i<max_skmer_size ; i++)
			// 	cout << skmer_buffer[i] << ", ";
			// cout << "]" << endl;

			if (position + 1 >= m)
			{
				// New minimizer
				if (candidate < minimizer)
				{
					// Output the previous kmer
					if (position >= k)
						saver->save_skmer(minimizer, skmer_buffer, max_skmer_size, skmer_start, position-1);

					// Save everything for new minimizer
					minimizer = candidate;
					minimizer_pos = position;
					if (position + 1 >= k)
						skmer_start = position - k + 1;
					// cout << "new\t" << candidate << '\t' << uint2kmer(candidate, m) << endl;
				}
				// Outdated minimizer
				else if (position - minimizer_pos > k - m) {
					// Output the previous skmer
					saver->save_skmer(minimizer, skmer_buffer, max_skmer_size, skmer_start, position-1);

					// Reinit the minimizer
					minimizer = 1ul << (2 * m);
					skmer_start = position - k + 1;

					// Search in near past mmers for the new minimizer
					for (size_t candidate_pos=position-k+m ; candidate_pos<=position ; candidate_pos++)
					{
						if (minimizer > candidates[candidate_pos % nb_mmers_in_kmer])
						{
							minimizer = candidates[candidate_pos % nb_mmers_in_kmer];
							minimizer_pos = candidate_pos;
						}
					}
					// cout << "outdated\t" << minimizer << '\t' << uint2kmer(minimizer, m) << endl;
				}
			}

			candidates[position % nb_mmers_in_kmer] = candidate;
			skmer_buffer[position % max_skmer_size] = c;
			position += 1;
		}
		c = fasta.next_char();
	}
	cout << "--------------- enumerating over -----------------" << endl;

	// Save the last superkmer
	saver->save_skmer(minimizer, skmer_buffer, max_skmer_size, skmer_start, position-1);
	delete saver;

	delete[] candidates;
	delete[] skmer_buffer;
}


int main(int argc, char const *argv[])
{
	compute_skmers(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));

	return 0;
}