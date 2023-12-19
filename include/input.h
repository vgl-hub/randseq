struct UserInputRandseq : UserInput {
    
    std::vector<std::string> outFiles; // output files
    float gcContent = 0.5, mutationRate = 0.0, readCoverage = 30.0;
    uint64_t sequenceLength, avgReadLen = 150;
    uint32_t randSeed = 1;

};

class Input {
    
    UserInputRandseq userInput;
    std::shared_ptr<std::istream> stream;
    
    char bases[4] = {'A','C','G','T'};
    float cumFreq[4] = {0.0, 0.0, 0.0, 0.0};
    
    
public:
    
    void load(UserInputRandseq userInput);

    void execute();
    
    int findCeil(float r);
    
};
