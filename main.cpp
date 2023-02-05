#include <string>
#include <fstream>
#include <iostream>

using namespace std;


class FastaStream
{
public:
	ifstream fs;

	FastaStream (string filename) {
		this->fs.open(filename);
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
		if (! this->has_next())
			return '\0';

		char c; bool done = false;
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
			while (c != '\n')
				this->fs.get(c);
			return '\0';
		default:
			return this->next_char();
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


#include <sstream>
#include <algorithm>
string uint2kmer (uint64_t val, size_t k) {
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


void compute_skmers (string filename, string outdir, size_t k, size_t m) {
	FastaStream fasta(filename);
	KmerStream stream(m);
	// Number of nucleotides to read before the first kmer
	size_t position = 0;
	uint64_t minimizer = 1ul << (2 * m);
	size_t minimizer_pos = 0;
	const size_t nb_mmers_in_kmer = k - m + 1;
	uint64_t * candidates = new uint64_t[nb_mmers_in_kmer];

	while (fasta.has_next())
	{
		char c = fasta.next_char();
		// New sequence
		if (c == '\0')
		{
			position = 0;
			minimizer = 1ul << (2 * m);
		}
		// New char in current sequence
		else
		{
			uint64_t candidate = stream.next_kmer(c);
			candidates[position % nb_mmers_in_kmer] = candidate;
			if (position + 1 >= m)
			{
				// New minimizer
				if (candidate < minimizer)
				{
					minimizer = candidate;
					minimizer_pos = position;
					cout << "new\t" << candidate << '\t' << uint2kmer(candidate, m) << endl;
				}
				// Outdated minimizer
				else if (position - minimizer_pos > k - m) {
					minimizer = 1ul << (2 * m);
					for (size_t candidate_pos=position-k+m ; candidate_pos<=position ; candidate_pos++)
					{
						if (minimizer > candidates[candidate_pos % nb_mmers_in_kmer])
						{
							minimizer = candidates[candidate_pos % nb_mmers_in_kmer];
							minimizer_pos = candidate_pos;
						}
					}
					cout << "outdated\t" << minimizer << '\t' << uint2kmer(minimizer, m) << endl;
				}
			}

			position += 1;
		}
	}

	delete[] candidates;
}


int main(int argc, char const *argv[])
{
	compute_skmers("data/small.fa", "data/small_out/", 5, 3);

	return 0;
}