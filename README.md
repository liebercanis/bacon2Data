code to analyze run2 bacon data

on branch RunTwo:

> git clone https://github.com/liebercanis/bacon2Data.git
> cd bacon2Data/

> git checkout -b runTwo
    Switched to a new branch 'runTwo'
>git branch --set-upstream-to=origin/runTwo runTwo
> git pull
cd to bobj
>rm *.o;
> make clean; make
cd to compiled
>make clean; make

main routine anaCRun.cc called by anac1.cc (main) requires environment variable BOBJ

then do, for example
>root "post.C(\"11_26_2023\")"
