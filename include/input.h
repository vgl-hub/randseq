#ifndef INPUT_H
#define INPUT_H

struct UserInputRandseq : UserInput {
    
    std::vector<std::string> outFiles; // output files
    float gcContent = 0.5, mutationRate = 0.0, readCoverage = 0;
    uint64_t sequenceLength, avgReadLen = 150;
    uint32_t randSeed = 1;

};

class Input {
    
    UserInputRandseq userInput;
    std::shared_ptr<std::istream> stream;
    StreamObj streamObj;
    
    char bases[4] = {'A','C','G','T'};
    float cumFreq[4] = {0.0, 0.0, 0.0, 0.0};
    
    std::vector<std::string> referenceCorrect;
    std::vector<std::pair<std::string, uint64_t>> errorPositions;
    
    
public:
    
    void load(UserInputRandseq userInput);
    
    void loadGenome(InSequences& inSequences);

    void execute();
    
    int findCeil(float r);
    
    void report();
    
    void writeFasta();
    
    void writeFastq();
    
};

#endif /* INPUT_H */
