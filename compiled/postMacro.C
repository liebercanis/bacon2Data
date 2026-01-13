
/*
look at  postAna root file
*/

// pass bit failures hex
enum FAILURECODES
{
    PASS = 0,
    BASEFAIL = 0x1,
    EARLYCUT = 0x2,
    FIRSTTIME = 0x4,
    COSMIC = 0x8,
    GAMMA = 0x10,
    TRIGFAIL = 0x20,
    TRIANGLE = 0x40, // 2^6
    TOTALCODES = 2 * TRIANGLE
};

enum
{
    FAILBITS = 8
};

std::string sdate;
TFile *fin;
TFile *fout;
std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number

string currentDate()
{
    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    char output[30];
    strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
    return string(output);
}

// collect trigger bits by file
void getTriggerBits()
{
    bitCount.resize(TOTALCODES / 2);

    TIter next(fin->GetListOfKeys());
    TKey *key;
    int ifile = 0;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());

        if (!cl->InheritsFrom("TH1D"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();

        if (!TString(h->GetName()).Contains("EventPassFile"))
            continue;

        // cout << h->GetName() << endl;

        for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
            bitCount[ibin].push_back(h->GetBinContent(ibin + 1));

        bitFile.push_back(++ifile);
    }

    cout << " #bits   " << bitCount.size() << "  # files " << bitCount[0].size() << endl;

    // normalize to number of files
    for (unsigned ibit = 0; ibit < bitCount[0].size(); ++ibit)
    {
        for (unsigned ifile = 0; ifile < bitCount[ibit].size(); ++ifile)
            bitCount[ibit][ifile] = bitCount[ibit][ifile] / double(bitFile.size());
    }

    // make graph for this bit
    for (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
    {
        int sum = std::accumulate(bitCount[ibit].begin(), bitCount[ibit].end(), 0.0);
        printf(" sum bit %i # name %s  sum %i  \n", ibit, codeNames[ibit].Data(), sum);
        if (sum > 1000)
        {
            printf(" graph for bit %i # name %s files %lu sum %i  \n", ibit, codeNames[ibit].Data(), bitCount[ibit].size(), sum);
            TGraph *gBit = new TGraph(bitFile.size(), &bitFile[0], &bitCount[ibit][0]);
            gBit->SetName(Form("TrigFailuresBit%i%s", ibit, codeNames[ibit].Data()));
            gBit->SetTitle(Form("TrigFailures for bit %i %s per file ", ibit, codeNames[ibit].Data()));
            gBit->SetMarkerStyle(21);
            gBit->GetYaxis()->SetTitle("falures per file ");
            gBit->GetXaxis()->SetTitle("file number ");
            fout->Append(gBit);
        }
    }
}

void postMacro()
{
    failCodes.resize(FAILBITS);
    failCodes[0] = PASS;
    failCodes[1] = BASEFAIL;
    failCodes[2] = EARLYCUT;
    failCodes[3] = FIRSTTIME;
    failCodes[4] = COSMIC;
    failCodes[5] = GAMMA;
    failCodes[6] = TRIGFAIL;
    failCodes[7] = TRIANGLE;
    bitNames.resize(FAILBITS);
    bitNames[0] = TString("pass");
    bitNames[1] = TString("baseline");
    bitNames[2] = TString("earlycut");
    bitNames[3] = TString("firsttime");
    bitNames[4] = TString("cosmic");
    bitNames[5] = TString("gamma");
    bitNames[6] = TString("trigger");
    bitNames[7] = TString("triangle");

    for (unsigned ic = 0; ic < TOTALCODES; ++ic)
        codeNames.push_back(TString("mixed"));
    codeNames[PASS] = TString("pass");
    codeNames[BASEFAIL] = TString("baseline");
    codeNames[TRIANGLE] = TString("triangle");
    codeNames[EARLYCUT] = TString("earlycut");
    codeNames[FIRSTTIME] = TString("firsttime");
    codeNames[COSMIC] = TString("cosmic");
    codeNames[GAMMA] = TString("gamma");

    // append to mixed names
    for (int ic = 0; ic < TOTALCODES; ++ic)
    {
        for (int ibit = 0; ibit < FAILBITS; ++ibit)
            if (ic & failCodes[ibit])
                codeNames[ic] += bitNames[ibit];
    }

    sdate = currentDate();
    printf(" making cans on %s \n", sdate.c_str());

    // put in explicit file name and get tag
    TString fileName("post-10_06_2025-10_06_2025-969976.root");
    TString tag = TString(fileName(fileName.First("-") + 1, 21));
    cout << " gains from file " << fileName << " with date tag " << tag << endl;

    // open file
    fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
        return;
    fout = new TFile("postMacro.root", "recreate");

    printf("open file %s date %s \n", fileName.Data(), sdate.c_str());
    getTriggerBits();
    fout->Write();
}