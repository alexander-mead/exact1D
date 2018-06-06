I started trying to add library support to all of this shit. I sort of half did it and then gave up. I did not test anything. This was probably very stupid.

At the moment (05/06/2018)

- I think preIC1D is probably unnecessary and could be deleted
- The RNG does not seem to work when real-8 or real-16 precision is specified. Exact1D needs very high precision to work
- I wanted to make the data from exact1D and PM1D compatible, but I have not done this
- I made a compile script ./compile.sh $1 where $1 is the name of the code to compile (preIC, IC1D, PM1D, exact1D)