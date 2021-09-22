#ifndef PROCESSLIST_H
#define PROCESSLIST_H

#include "Process.h"

#include <TString.h>
#include <TFile.h>
#include <THStack.h>

class ProcessList {
    private:
        Process* head = nullptr;
        Process* tail = nullptr;

    public:
        ProcessList() = default;
        ~ProcessList();

        void addProcess(TString& name, int color, TFile* inputfile);

        Process* getHead() {return head;}
        Process* getTail() {return tail;}

        void fillStack(THStack* stack, TString& histogramID);

};

#endif