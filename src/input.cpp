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

#include "stream-obj.h"
#include "fastx.h"
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

    lg.verbose("Sequence length: " + std::to_string(userInput.sequenceLength));
    lg.verbose("GC content: " + std::to_string(userInput.gcContent));
    lg.verbose("Average read length: " + std::to_string(userInput.avgReadLen));
    lg.verbose("Mutation rate: " + std::to_string(userInput.mutationRate));
    lg.verbose("Read coverage: " + std::to_string(userInput.readCoverage));
    lg.verbose("Random seed: " + std::to_string(userInput.randSeed));
    lg.verbose("Cumulative frequency [A,C,G,T]: " + std::to_string(cumFreq[0]) + "," + std::to_string(cumFreq[1]) + "," + std::to_string(cumFreq[2]) + "," + std::to_string(cumFreq[3]));
    
    srand(userInput.randSeed);
    report();
    
    //consolidate log
//    for (auto it = logs.begin(); it != logs.end(); it++) {
//     
//        it->print();
//        logs.erase(it--);
//        if(verbose_flag) {std::cerr<<"\n";};
//        
//    }
    
}

void Input::loadGenome(InSequences& inSequences) {
    
    if (userInput.inSequence.empty()) {return;}
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    stream = streamObj.openStream(userInput, 'f'); // open file
    
    if (stream) {
        
        switch (stream->peek()) {
                
            case '>': {
                
                stream->get();
                
                while (getline(*stream, newLine)) {
                    
                    h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    seqHeader = h;
                    
                    if (c != NULL)
                        seqComment = std::string(c);
                    
                    std::string* inSequence = new std::string;
                    
                    getline(*stream, *inSequence, '>');
                    
                    lg.verbose("Individual fasta sequence read.");
                    
                    Sequence* sequence = new Sequence {seqHeader, seqComment, inSequence};
                    
                    if (sequence != NULL) {
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence, -1);
                        
                        seqPos++;
                        
                    }
                    
                }
                
                break;
            }
            case '@': {
                
                while (getline(*stream, newLine)) { // file input
                    
                    newLine.erase(0, 1);
                    
                    h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    seqHeader = h;
                    
                    if (c != NULL) {
                        
                        seqComment = std::string(c);
                        
                    }
                    
                    std::string* inSequence = new std::string;
                    getline(*stream, *inSequence);
                    
                    getline(*stream, newLine);
                    
                    std::string* inSequenceQuality = new std::string;
                    getline(*stream, *inSequenceQuality);
                    
                    Sequence* sequence = new Sequence {seqHeader, seqComment, inSequence, inSequenceQuality};
                    
                    if (sequence != NULL) {
                        
                        sequence->seqPos = seqPos; // remember the order
                    
                        inSequences.appendSequence(sequence, -1);
                        
                        seqPos++;
                        
                    }
                    
                }
                
                break;
                
            }
            default: {
                
                readGFA(inSequences, userInput, stream);
                
            }
            
        }
        
        lg.verbose("End of file.");
            
    }else{

        fprintf(stderr, "Stream not successful: %s.", userInput.inSequence.c_str());
        exit(1);

    }

    jobWait(threadPool);
    
    inSequences.updateStats(); // compute summary statistics

}

void Input::report() { // generates the output from the program
    
    const static phmap::flat_hash_map<std::string,int> string_to_case { // different outputs available
        {"fasta",1},
        {"bed",2},
        {"fastq",3}
    };
    
    std::string ext = "stdout";
    
    //if (userInput.outFile != "")
    //    ext = getFileExt("." + userInput.outFile);
    
    //lg.verbose("Writing ouput: " + userInput.outFile);
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
        
        default:
        case 1: { // .fasta
            
            writeFasta();
            if (userInput.readCoverage != 0)
                writeFastq();
            
            break;
            
        }
            
        case 2: { // .fastq
            
            //writeFastq();
            
            break;
            
        }
            
        case 3: { // .bed
            
            //writeBed();
            
            break;
            
        }
            
    }
    
}

void Input::writeFasta() {
    
    float r;
    char base, errorBase;
    std::string referenceCorrectSeq;
    
    std::ofstream referenceCorrectFile;
    referenceCorrectFile.open("referenceCorrect.fasta");
    
    std::ofstream referenceErrorFile;
    referenceErrorFile.open("referenceError.fasta");
    
    std::ofstream errorVcfFile;
    errorVcfFile.open("errors.vcf");
    
    std::string header = "Reference";
    
    if (userInput.inSequence.empty()) {
        
        referenceCorrectFile << ">Reference\n";
        referenceErrorFile << ">Reference with errors\n";
        
        for (uint64_t i = 0; i < userInput.sequenceLength; ++i) { // generate reference
            
            r = newRand();
            base = bases[findCeil(r)];
            referenceCorrectFile<<base;
            referenceCorrectSeq.push_back(base);
            
            if (newRand() <= userInput.mutationRate) {
                
                errorBase = bases[findCeil(newRand())];
                
                referenceErrorFile<<errorBase;
                errorVcfFile<<header<<"\t"<<i<<"\t"<<base<<"\t"<<errorBase<<std::endl;;
                
            }else{
                
                referenceErrorFile<<base;
                
            }
            
        }
        
        referenceCorrectFile<<std::endl;
        referenceErrorFile<<std::endl;
        referenceCorrect.push_back(referenceCorrectSeq);
        
    }else{
        
        InSequences genome; // initialize sequence collection object
            
        if (!userInput.inSequence.empty()) {
            lg.verbose("Loading input sequences");
            loadGenome(genome); // read input genome
            lg.verbose("Sequences loaded");
        }
    
        genome.sortPathsByOriginal();
        
        std::vector<InPath> inPaths = genome.getInPaths();
        std::vector<InSegment*> *inSegments = genome.getInSegments();
        std::vector<InGap> *inGaps = genome.getInGaps();
        
        for (InPath& path : inPaths) {
            
            unsigned int cUId = 0, gapLen = 0;
            
            std::vector<PathComponent> pathComponents = path.getComponents();
            
            uint64_t absPos = 0;
            
            referenceCorrectFile<<">"<<path.getHeader()<<std::endl;
            referenceErrorFile<<">"<<path.getHeader()<<std::endl;
            
            header = path.getHeader();
            
            referenceCorrectSeq.clear();
            
            for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
                
                cUId = component->id;
                
                if (component->type == SEGMENT) {
                    
                    auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                    
                    std::string sequence = (*inSegment)->getInSequence();
                    
                    if (component->orientation == '+') {
                        
                        for (uint64_t i = 0; i < (*inSegment)->getSegmentLen(); ++i) {
                            
                            base = sequence[i];
                            referenceCorrectFile<<base;
                            referenceCorrectSeq.push_back(base);
                            
                            r = newRand();
                            
                            if (newRand() <= userInput.mutationRate) {
                                
                                errorBase = bases[findCeil(newRand())];
                                
                                referenceErrorFile<<errorBase;
                                errorVcfFile<<header<<"\t"<<absPos<<"\t"<<base<<"\t"<<errorBase<<std::endl;;
                                
                            }else{
                                referenceErrorFile<<base;
                            }
                            ++absPos;
                        }
                    }else{
                        // GFA not handled yet
                    }
                }else if (component->type == GAP){
                    
                    auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                    
                    gapLen += inGap->getDist(component->start - component->end);
                    absPos += gapLen;
                    
                }else{} // need to handle edges, cigars etc
            }
            referenceCorrectFile<<std::endl;
            referenceErrorFile<<std::endl;
            referenceCorrect.push_back(referenceCorrectSeq);
        }
    }
    referenceCorrectFile.close();
    referenceErrorFile.close();
    errorVcfFile.close();
}

void Input::writeFastq() {
    
    uint16_t pass = round(userInput.avgReadLen / userInput.readCoverage);
    
    std::ofstream readsFile;
    readsFile.open ("reads.fastq");
    
    for (std::string seq : referenceCorrect) {
        
        for (uint64_t i = 0; i < seq.size(); ++i) {
            
            if ((i % pass == 0) && (i + userInput.avgReadLen < seq.size())) {
                
                std::cout<<pass<<std::endl;
                
                readsFile<< "@read " << i << "\n"
                            << seq.substr(i, userInput.avgReadLen) << "\n"
                            << "+\n"
                            << std::string(userInput.avgReadLen, '!') << "\n";
                
            }
            
        }
        
    }

    readsFile.close();
    
}
