Instructions to run mcfm (note to self).

Install mcfm, see the [mcfm home page](http://mcfm.fnal.gov/)

Then

```
cd Bin
# edit input.DAT: set proc=261, beam energy, pp, etc.
./mcfm
h2root file.rz file.root
```

See the mcfm manual for
1. the processes that you should include
2. the way in which the particles are labeled
