
/* get rawBr */
unsigned trigStart = 600;
unsigned trigEnd = 1200;
TFile *fout;
TFile *fin;
TTree *rawTree;
int passBit;
// vectors for gains
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
//
std::map<int, int> chanMap;
vector<TBRawEvent *> rawBr;
TBEventData *eventData;
TBEventData *rawEventData;
TBRun *tbrun;
TNtuple *ntChan;
TNtuple *ntChanSum;
vector<TH1D *> baseHist;
vector<TH1D *> noiseHist;
vector<TH1D *> skewHist;
vector<TH1D *> sumWave;
vector<TH1D *> sumHitWave;
vector<TH1D *> sumPeakWave;
vector<TH1D *> valHist;
vector<TH1D *> sumWaveB;
vector<TH1D *> valHistB;

vector<TH1D *> hQSum;
vector<TH1D *> hQPeak;
vector<TH1D *> hQPEShape;
TH1D *hEvBaseWave;
vector<TH1D *> hEvGaus;
vector<TH1D *> hChannelGaus;
TH1D *evCount;
TH1D *histQSum;
TH1D *hEventPass;
TH1D *histHitCount;
TH1D *hSumPMT;
TH1D *threshHist;
TH2D *threshValueHist;
TH1D *crossHist;
TH1D *cosmicCut1;
TH1D *cosmicCut2;
// TH1D *histQPE;
TH1D *histQPrompt;
vector<double> channelSigmaValue;
vector<double> channelSigma;
vector<double> channelSigmaErr;
vector<double> digi;
vector<double> ddigi;
vector<double> hdigi;
std::vector<unsigned> thresholds;
std::vector<unsigned> crossings;
std::vector<unsigned> crossingBin;
std::vector<double> crossingTime;
vector<double> slope;
vector<double> eslope;
vector<double> chan;
vector<double> echan;

unsigned getBranches()
{
    TObjArray *brList = rawTree->GetListOfBranches();
    brList->ls();
    TString cname;
    TIter next(brList);
    TBranch *aBranch = NULL;
    while ((aBranch = (TBranch *)next()))
    {
        TString s(aBranch->GetName());
        if (s != TString("eventData"))
        {
            int ichan = TString(s(s.Last('n') + 1, s.Length())).Atoi();
            // rawTree->GetBranch(aBranch->GetName())->SetAutoDelete(kTRUE);
            cout << s << "  " << aBranch->GetName() << " return val =  " << rawTree->SetBranchAddress(aBranch->GetName(), &rawBr[ichan]) << endl;
        }
    }
    return rawBr.size();
}

bool openFile(TString theFile)
{
    // open input file and make some histograms
    TString fileName;
    fileName.Form("rootData/%s", theFile.Data());
    printf(" looking for file %s\n", fileName.Data());

    bool exists = false;
    FILE *aFile;
    aFile = fopen(fileName.Data(), "r");
    if (aFile)
    {
        fclose(aFile);
        exists = true;
    }
    if (!exists)
    {
        printf(" couldnt open file %s\n", fileName.Data());
        return false;
    }

    fin = new TFile(fileName, "readonly");
    printf(" opened file %s\n", fileName.Data());
    rawTree = NULL;
    fin->GetObject("RawTree", rawTree);
    if (!rawTree)
    {
        printf(" no RawTree in file %s\n", fileName.Data());
        return false;
    }
    cout << "  RawTree has " << rawTree->GetEntries() << " entries " << endl;
    rawEventData = new TBEventData();
    rawTree->SetBranchAddress("eventData", &rawEventData);
    if (!rawEventData)
    {
        printf(" eventData not found in file  %s\n", fileName.Data());
        return false;
    }
    printf(" rawTree has %u channels stored in rawBr \n", getBranches());
    for (unsigned i = 0; i < rawBr.size(); ++i)
        printf(" branch %s chan %i \n", rawBr[i]->GetName(), i);

    return true;
}

void anaTest()
{
    rawBr.clear();

    for (int ichan = 0; ichan < 13; ++ichan)
    {
        TBRawEvent *rawEv = new TBRawEvent(ichan);
        rawEv->rdigi.resize(7500);
        rawEv->SetName(Form("rawChan%i", ichan));
        rawBr.push_back(rawEv);
    }
    openFile("run-11_26_2023-file.root");
    //cout << "read gains from file " << gainFileName << endl;
    //readGains(gainFileName);

    // need to fill rawBr[0]->rdigi.size()
    printf("Read zeroth entry from tree \n");
    if (!rawTree)
    {
        printf("EEEEEE rawTree is null!!!!!\n");
        return 0;
    }
    cout << " RawTree still has has " << rawTree->GetEntries() << " entries " << endl;
    cout << " rawTree return " << rawTree->GetEntry(0) << endl;
    printf("got rawTree entry 0 \n");
    printf("\n\n\t\t >>>>>>>>> start of file %i %i %i : %i <<<<<<<<<<<< \n", rawEventData->day, rawEventData->mon, rawEventData->year, rawEventData->hour);
}