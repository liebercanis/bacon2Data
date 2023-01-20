TFile *fin;
vector<TTree *> treeList;
vector<int> chanList;
vector<TBRawEvent *> rawBr;
vector<TH1D *> sumWave;
TCanvas *can;

unsigned getTrees()
{
    treeList.clear();
    chanList.clear();
    //
    TIter next(fin->GetListOfKeys());
    TKey *key;
    unsigned count = 0;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TTree"))
            continue;
        TTree *tree = (TTree *)key->ReadObj();
        // cout << tree->GetName() << " " << key->GetCycle() << endl;
        if (tree->GetEntries() < 1)
            continue;
        // pick up only last in cycle
        int ichan = TString(TString(tree->GetName())(5, 2)).Atoi();
        if (find(chanList.begin(), chanList.end(), ichan) != chanList.end())
            continue;
        cout << count << " tree name " << tree->GetName() << " chan " << ichan << endl;
        treeList.push_back(tree);
        chanList.push_back(ichan);
        ++count;
    }
    return treeList.size();
}

/* get rawBr */
void getBranches()
{
    for (unsigned it = 0; it < treeList.size(); ++it)
    {
        rawBr[it] = NULL;
        TObjArray *brList = treeList[it]->GetListOfBranches();
        TIter next(brList);
        TBranchElement *aBranch = NULL;
        while ((aBranch = (TBranchElement *)next()))
        {
            // cout << "tree " << it << " " << aBranch->GetName() << endl;
            TBRawEvent *chan = (TBRawEvent *)aBranch->GetObject();
            if (!chan)
            {
                cout << " \t\t !!!! NULL " << it << endl;
                continue;
            }
            rawBr[it] = chan;
            cout << it << " rawBr name  ==>  " << rawBr[it]->GetName() << endl;
        }
    }
}

void getEvent(Long64_t entry)
{
    for (unsigned it = 0; it < treeList.size(); ++it)
        treeList[it]->GetEntry(entry);
    getBranches();
}


void plots()
{
    TIter next(fin->GetListOfKeys());
    TKey *key;
    TCanvas *can;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1"))
            continue;
        TH1I *h = (TH1I *)key->ReadObj();
        TString hname(h->GetName());
        cout << hname << endl;
        TString canName;
        canName.Form("plot-%s", hname.Data());
        can = new TCanvas(canName, canName);
        can->SetLogy();
        h->Draw();
        can->Print(".pdf");
    }
}

void plotWave(unsigned ichan, Long64_t entry)
{
    getEvent(entry);
    // find corresponding tree
    unsigned ib = 1000;
    for (unsigned it = 0; it < chanList.size(); ++it)
        if (chanList[it] == ichan)
            ib = it;
    if(ib==1000){
        printf("no tree for chan %u \n",ichan);
        return;
    }
    printf("tree for chan %u numbered %u \n",ichan,ib);

    unsigned baseLength = 30;
    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
        base += (double)rawBr[ib]->rdigi[j];
    }
    base /= double(baseLength);
    sumWave[ib]->Reset("ICES");

    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
        double val = (double)rawBr[ib]->rdigi[j] - base;

        sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + val);
    }

    TString canName;
    canName.Form("event%lldchan%u",entry,ib);
    can = new TCanvas(canName,canName);
    sumWave[ib]->SetTitle(canName);
    sumWave[ib]->Draw();
}
/* main */
void plotSisEvent(TString file = "data/rootData/run-01_12_2023-nev0.root")
{
    fin = new TFile(file, "readonly");
    fin->ls();
    getTrees();
    cout <<" size of chanList " << chanList.size() << endl;
    rawBr.resize(chanList.size());
    getEvent(0);

    printf("make hist \n");

    TFile *fout = new TFile("plotSisEvent.root", "recreate");

    for (unsigned i = 0; i < rawBr.size(); ++i)
    {
        unsigned ichan = chanList[i];
        sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
    }
}