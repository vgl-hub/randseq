#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <getopt.h>
#include <vector>

#include "log.h"
#include "uid-generator.h"
#include "bed.h"
#include "global.h"
#include "struct.h"
#include "functions.h"
#include "threadpool.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"
#include "input.h"
#include "main.h"

std::string version = "0.0.1";

std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int verbose_flag;
int tabular_flag;
int maxThreads = 0;
int cmd_flag;
std::vector<Log> logs;
UserInputRandseq userInput;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;

std::string getArgs(char* optarg, unsigned int argc, char **argv) {
    
    std::string cmd;
    bool record = false;

    for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
        
        if (optarg != argv[arg_counter] && !record) {
            continue;
        }else{
            record = true;
            if(optarg != argv[arg_counter]){
                cmd += ' ';
                cmd += argv[arg_counter];
            }
        }
    }
    
    return cmd;
    
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
//    std::string cmd;
    
    bool isPipe = false; // to check if input is from pipe
    
    if (argc == 1) { // mytool with no arguments
            
        printf("randseq [command]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"average-read-length", required_argument, 0, 'a'},
        {"input-sequence", required_argument, 0, 'f'},
        {"gc-content", required_argument, 0, 'g'},
        {"sequence-length", required_argument, 0, 'l'},
        {"mutation-rate", required_argument, 0, 'm'},
        {"read-coverage", required_argument, 0, 'r'},
        {"rand-seed", required_argument, 0, 's'},
        
        {"threads", required_argument, 0, 'j'},
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    const static std::unordered_map<std::string,int> tools{
        {"something",1},
        {"somethingelse",2}
    };
    
    while (arguments) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:a:f:g:l:m:r:s:j:vh",
                        long_options, &option_index);
        
        if (c == -1) { // exit the loop if run out of options
            break;
            
        }

        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'b':
                        break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default: // handle positional arguments
                                
                switch (tools.count(optarg) ? tools.at(optarg) : 0) {
                    case 1:
                        //cmd = "gfastats/build/bin/mytool" + getArgs(optarg, argc, argv);;
                        
                        arguments = false;
                        
                        break;
                    case 2:
                        //cmd = "gfalign/build/bin/mytool2" + getArgs(optarg, argc, argv);;
                        
                        arguments = false;
                        
                        break;
                        
                }
                
            case 0: // case for long options without short options
                
//                if (strcmp(long_options[option_index].name,"line-length") == 0)
//                  splitLength = atoi(optarg);
                
                break;
                
            case 'a': // average read length
                userInput.avgReadLen = atoi(optarg);
                break;

            case 'f': // input sequence
                
                if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    userInput.pipeType = 'f'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    userInput.inSequence = optarg;
                    userInput.stats_flag = 1;
                    
                }
                    
                break;
                
            case 'g': // gc content
                userInput.gcContent = atof(optarg);
                break;
                
            case 'l': // sequence length
                userInput.sequenceLength = atoi(optarg);
                break;
                
            case 'm': // mutation-rate
                userInput.mutationRate = atof(optarg);
                break;
                
            case 'r': // read coverage
                userInput.readCoverage = atof(optarg);
                break;
                
            case 's': // read coverage
                userInput.randSeed = atoi(optarg);
                break;
                
            case 'j': // max threads
                maxThreads = atoi(optarg);
                userInput.stats_flag = 1;
                break;
                
            case 'v': // software version
                printf("randseq v%s\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com\n");
                exit(0);
                
            case 'h': // help
                printf("randseq [command]\n");
                printf("\nOptions:\n");
                printf("\t-f --input-sequence sequence input file (fasta,gfa1/2).\n");
                printf("\t-a --average-read-length <int> (default: 150).\n");
                printf("\t-g --gc-content <float> (default: 0.5).\n");
                printf("\t-l --sequence-length <int>.\n");
                printf("\t-m --mutation-rate <float> (default: 0).\n");
                printf("\t-r --read-coverage <float> (default: 0).\n");
                printf("\t-s --rand-seed <int>.\n");
                printf("\t-v --version software version.\n");
                printf("\t--verbose verbose output.\n");
                printf("\t--cmd print $0 to stdout.\n");
                exit(0);
        }
        
        if    (argc == 2 || // handle various cases in which the output should include summary stats
              (argc == 3 && pos_op == 2) ||
              (argc == 4 && pos_op == 3)) {
            
        }
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    Input in;
    
    threadPool.init(maxThreads); // initialize threadpool
    
    in.load(userInput); // load user input
    
    lg.verbose("Loaded user input");
    
    in.execute(); // execute tool
    
    threadPool.join(); // join threads
    
//    std::cout<<"Invoking: "<<cmd<<std::endl;
//    std::system(cmd.c_str());
    
    exit(EXIT_SUCCESS);
    
}
