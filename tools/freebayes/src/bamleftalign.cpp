#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>

#include "Fasta.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
//#include "IndelAllele.h"
#include "LeftAlign.h"

#ifdef VERBOSE_DEBUG
#define DEBUG(msg) \
    if (debug) { cerr << msg; }
#else
#define DEBUG(msg)
#endif

using namespace std;
using namespace BamTools;

void printUsage(char** argv) {
    cerr << "usage: [BAM data stream] | " << argv[0] << " [options]" << endl
         << endl
         << "Left-aligns and merges the insertions and deletions in all alignments in stdin." << endl
         << "Iterates until each alignment is stable through a left-realignment step." << endl
         << endl
         << "arguments:" << endl
         << "      -f --fasta-reference FILE   FASTA reference file to use for realignment (required)" << endl
         << "      -d --debug             Print debugging information about realignment process" << endl
         << "      -s --suppress-output   Don't write BAM output stream (for debugging)" << endl
         << "      -m --max-iterations N  Iterate the left-realignment no more than this many times" << endl
         << "      -c --compressed        Write compressed BAM on stdout, default is uncompressed" << endl;
}

int main(int argc, char** argv) {

    int c;

    FastaReference reference;
    bool has_ref = false;
    bool suppress_output = false;
    bool debug = false;
    bool isuncompressed = true;

    int maxiterations = 50;
    
    if (argc < 2) {
        printUsage(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"debug", no_argument, 0, 'd'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"max-iterations", required_argument, 0, 'm'},
            {"suppress-output", no_argument, 0, 's'},
            {"compressed", no_argument, 0, 'c'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hdcsf:m:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c) {

            case 'f':
                reference.open(optarg); // will exit on open failure
                has_ref = true;
                break;
     
            case 'm':
                maxiterations = atoi(optarg);
                break;

            case 'd':
                debug = true;
                break;

            case 's':
                suppress_output = true;
                break;

            case 'c':
                isuncompressed = false;
                break;

            case 'h':
                printUsage(argv);
                exit(0);
                break;
              
            case '?':
                printUsage(argv);
                exit(1);
                break;
     
              default:
                abort();
                break;
        }
    }

    if (!has_ref) {
        cerr << "no FASTA reference provided, cannot realign" << endl;
        exit(1);
    }

    BamReader reader;
    if (!reader.Open("stdin")) {
        cerr << "could not open stdin for reading" << endl;
        exit(1);
    }

    BamWriter writer;

    if (isuncompressed) {
        writer.SetCompressionMode(BamWriter::Uncompressed);
    }

    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    BamAlignment alignment;

    while (reader.GetNextAlignment(alignment)) {

            DEBUG("---------------------------   read    --------------------------" << endl);
            DEBUG("| " << referenceIDToName[alignment.RefID] << ":" << alignment.Position << endl);
            DEBUG("| " << alignment.Name << ":" << alignment.GetEndPosition() << endl);
            DEBUG("| " << alignment.Name << ":" << (alignment.IsMapped() ? " mapped" : " unmapped") << endl);
            DEBUG("| " << alignment.Name << ":" << " cigar data size: " << alignment.CigarData.size() << endl);
            DEBUG("--------------------------- realigned --------------------------" << endl);

            // skip unmapped alignments, as they cannot be left-realigned without CIGAR data
            if (alignment.IsMapped()) {

                int endpos = alignment.GetEndPosition();
                int length = endpos - alignment.Position + 1;
                if (alignment.Position >= 0 && length > 0) {
                    if (!stablyLeftAlign(alignment,
                                reference.getSubSequence(
                                    referenceIDToName[alignment.RefID],
                                    alignment.Position,
                                    length),
                                maxiterations, debug)) {
                        cerr << "unstable realignment of " << alignment.Name
                             << " at " << referenceIDToName[alignment.RefID] << ":" << alignment.Position << endl
                             << alignment.AlignedBases << endl;
                    }
                }

            }

            DEBUG("----------------------------------------------------------------" << endl);
            DEBUG(endl);

        if (!suppress_output)
            writer.SaveAlignment(alignment);

    }

    reader.Close();
    if (!suppress_output)
        writer.Close();

    return 0;

}
