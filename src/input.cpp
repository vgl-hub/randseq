#include <stdlib.h>
#include <string>

#include <istream>
#include <fstream>
#include <sstream>

#include <parallel-hashmap/phmap.h>

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"

#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"
#include "input.h"

float newRand() {
    
    return (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
    
}

void Input::load(UserInputRandseq userInput) {
    
    this->userInput = userInput;
    cumFreq[0] = (1.0 - userInput.gcContent) / 2;
    cumFreq[1] = (userInput.gcContent) / 2 + cumFreq[0];
    cumFreq[2] = (userInput.gcContent) / 2 + cumFreq[1];
    cumFreq[3] = (1.0 - userInput.gcContent) / 2 + cumFreq[2];
    
}

int Input::findCeil(float r) {
    for (uint8_t i = 0; i < 4; ++i) {
        
        if (r <= cumFreq[i])
            return i;
        
    }
    return 4;
}
    
void Input::execute() {
    
    //if (userInput.inSequence.empty()) {return;}
    
    //threadPool.init(maxThreads); // initialize threadpool
    
    lg.verbose("Sequence length: " + std::to_string(userInput.sequenceLength));
    lg.verbose("GC content: " + std::to_string(userInput.gcContent));
    lg.verbose("Average read length: " + std::to_string(userInput.avgReadLen));
    lg.verbose("Mutation rate: " + std::to_string(userInput.mutationRate));
    lg.verbose("Read coverage: " + std::to_string(userInput.readCoverage));
    lg.verbose("Random seed: " + std::to_string(userInput.randSeed));
    lg.verbose("Cumulative frequency [A,C,G,T]: " + std::to_string(cumFreq[0]) + "," + std::to_string(cumFreq[1]) + "," + std::to_string(cumFreq[2]) + "," + std::to_string(cumFreq[3]));
    
    srand(userInput.randSeed);
    
    float r;
    char base;
    
    std::string referenceCorrect, referenceError;
    std::vector<std::string> reads;
    std::vector<uint64_t> errorPositions;
    
    uint16_t pass = round(userInput.avgReadLen / userInput.readCoverage);

    for (uint64_t i = 0; i < userInput.sequenceLength; ++i) { // generate reference
        
        r = newRand();
        base = bases[findCeil(r)];
        referenceCorrect.push_back(base);
        
        if (newRand() <= userInput.mutationRate) {
            
            referenceError.push_back(bases[findCeil(newRand())]);
            errorPositions.push_back(i);
            
        }else{
            
            referenceError.push_back(base);
            
        }
        
    }
    
    for (uint64_t i = 0; i < userInput.sequenceLength; ++i) {
        
        if ((i % pass == 0) && (i + userInput.avgReadLen < userInput.sequenceLength)) {
            
            reads.push_back(referenceCorrect.substr(i, userInput.avgReadLen));
            
        }
        
    }
    
    std::ofstream referenceCorrectFile;
    referenceCorrectFile.open ("referenceCorrect.fasta");
    referenceCorrectFile << ">Reference\n" << referenceCorrect << std::endl;
    referenceCorrectFile.close();
    
    std::ofstream referenceErrorFile;
    referenceErrorFile.open ("referenceError.fasta");
    referenceErrorFile << ">Reference with errors\n" << referenceError << std::endl;
    referenceErrorFile.close();
    
    std::ofstream readsFile;
    readsFile.open ("reads.fastq");
    for (uint64_t i = 0; i < reads.size(); ++i) {
        readsFile<< "@read " << i << "\n"
                    << reads[i] << "\n"
                    << "+\n"
                    << std::string(reads[i].size(), '!') << "\n";
    }
    readsFile.close();
    
    std::ofstream errorBedFile;
    errorBedFile.open ("errors.bed");
    for (uint64_t coord : errorPositions) {
        errorBedFile << "Reference\t" << coord << "\n";
    }
    errorBedFile.close();
    
    //jobWait(threadPool);
    
    //consolidate log
//    for (auto it = logs.begin(); it != logs.end(); it++) {
//     
//        it->print();
//        logs.erase(it--);
//        if(verbose_flag) {std::cerr<<"\n";};
//        
//    }
    
//    threadPool.join();
    
}
