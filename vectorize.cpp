#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <math.h>

using namespace std;

u_int64_t revComp(u_int64_t x, size_t sizeKmer)
{
    u_int64_t res = x;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2 * (32 - sizeKmer)));
}

unordered_map<int, double> kSizeMers(int size)
{
    unordered_map<int, double> kmers;

    for (u_int64_t i = 0; i < (int)pow(4, size); i++)
    {
        kmers[min(i, revComp(i, size))] = 0;
    }

    return kmers;
}

vector<double> getKmersProfile(string &line, int size, bool normalize, unordered_map<int, double> &KMERS)
{
    double total = 0;
    long len = 0;
    u_int64_t val = 0;
    unordered_map<int, double> kmers;
    kmers.insert(KMERS.begin(), KMERS.end());
    vector<double> stats(kmers.size(), 0);

    for (int i = 0; i < (int)line.length(); i++)
    {
        if (!(line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T'))
        {
            val = 0;
            len = 0;
            continue;
        }

        val = (val << 2);
        val = val & (int)(pow(2, 2 * size) - 1);
        val += (line[i] >> 1 & 3);
        len++;

        if (len == size)
        {
            // use val as the kmer for counting
            len--;
            kmers[min(val, revComp(val, size))]++;
            total++;
        }
    }

    int i = 0;
    for (unordered_map<int, double>::iterator it = kmers.begin(); it != kmers.end(); ++it)
    {
        if (normalize)
        {
            stats[i] = it->second / max(1.0, total);
        }
        else
        {
            stats[i] = it->second;
        }
        i++;
    }

    return stats;
}

void processBatch(vector<string> &lines, int kmerSize, int threads)
{
    vector<vector<double>> batchAnswers(lines.size());
    unordered_map<int, double> kmers = kSizeMers(kmerSize);

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < lines.size(); i++)
    {
        batchAnswers[i] = getKmersProfile(lines[i], kmerSize, true, kmers);
    }

    for (vector<double> batchAnswer : batchAnswers)
    {
        for (int j = 0; j < batchAnswer.size(); j++)
        {
            cout << batchAnswer[j];
            if (j < batchAnswer.size() - 1)
            {
                cout << "\t";
            }
        }
        cout << endl;
    }
}

int main(int argc, char **argv)
{
    string reads_path = argv[1];
    int kmerSize = stoi(argv[2]);
    int threads = stoi(argv[3]);
    int batchSize = stoi(argv[4]);

    ifstream reads(reads_path);
    vector<string> batch;
    string line;
    int lineNo = 0;

    while (getline(reads, line))
    {
        if (lineNo % 4 == 1)
        {
            batch.push_back(line);
        }

        lineNo++;

        if (batch.size() == batchSize)
        {
            processBatch(batch, kmerSize, threads);
            batch.clear();
        }
    }
    if (batch.size() > 0)
    {
        processBatch(batch, kmerSize, threads);
        batch.clear();
    }

    reads.close();

    return 0;
}
