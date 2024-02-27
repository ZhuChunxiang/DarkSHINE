// #ifndef rnd_h
// #define rnd_h

#include "t1.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class rnd : public t1 {
    public:
        rnd(TTree *tree=0):t1(tree){};
        ~rnd() {};
        virtual void Loop(string filename);
        void initialize(string filename);
        void execute();
        void finalize();
        void clear();

    private:
        TFile *m_file;
        vector<double>  *layer_Edep;
        vector<int>     *Photon_num;
        TBranch         *b_Photon_num;   
        TBranch         *b_layer_Edep;   
        TTree           *fChain;
        char            filename;

};
// #endif