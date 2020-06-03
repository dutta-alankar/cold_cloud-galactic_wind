# cold_cloud-galactic_wind
Script for hydrodynamic modelling cold clouds interacting with hot galactic winds using PLUTO Hydrodynamic solver (version 4.3)

```
 ▄████▄   ██▓     ▒█████   █    ██ ▓█████▄     ▄████▄   ██▀███   █    ██   ██████  ██░ ██  ██▓ ███▄    █   ▄████ 
▒██▀ ▀█  ▓██▒    ▒██▒  ██▒ ██  ▓██▒▒██▀ ██▌   ▒██▀ ▀█  ▓██ ▒ ██▒ ██  ▓██▒▒██    ▒ ▓██░ ██▒▓██▒ ██ ▀█   █  ██▒ ▀█▒
▒▓█    ▄ ▒██░    ▒██░  ██▒▓██  ▒██░░██   █▌   ▒▓█    ▄ ▓██ ░▄█ ▒▓██  ▒██░░ ▓██▄   ▒██▀▀██░▒██▒▓██  ▀█ ██▒▒██░▄▄▄░
▒▓▓▄ ▄██▒▒██░    ▒██   ██░▓▓█  ░██░░▓█▄   ▌   ▒▓▓▄ ▄██▒▒██▀▀█▄  ▓▓█  ░██░  ▒   ██▒░▓█ ░██ ░██░▓██▒  ▐▌██▒░▓█  ██▓
▒ ▓███▀ ░░██████▒░ ████▓▒░▒▒█████▓ ░▒████▓    ▒ ▓███▀ ░░██▓ ▒██▒▒▒█████▓ ▒██████▒▒░▓█▒░██▓░██░▒██░   ▓██░░▒▓███▀▒
░ ░▒ ▒  ░░ ▒░▓  ░░ ▒░▒░▒░ ░▒▓▒ ▒ ▒  ▒▒▓  ▒    ░ ░▒ ▒  ░░ ▒▓ ░▒▓░░▒▓▒ ▒ ▒ ▒ ▒▓▒ ▒ ░ ▒ ░░▒░▒░▓  ░ ▒░   ▒ ▒  ░▒   ▒ 
  ░  ▒   ░ ░ ▒  ░  ░ ▒ ▒░ ░░▒░ ░ ░  ░ ▒  ▒      ░  ▒     ░▒ ░ ▒░░░▒░ ░ ░ ░ ░▒  ░ ░ ▒ ░▒░ ░ ▒ ░░ ░░   ░ ▒░  ░   ░ 
░          ░ ░   ░ ░ ░ ▒   ░░░ ░ ░  ░ ░  ░    ░          ░░   ░  ░░░ ░ ░ ░  ░  ░   ░  ░░ ░ ▒ ░   ░   ░ ░ ░ ░   ░ 
░ ░          ░  ░    ░ ░     ░        ░       ░ ░         ░        ░           ░   ░  ░  ░ ░           ░       ░ 
░                                   ░         ░                                                                  
```

The work is in line with the method discussed in "The growth and entrainment of cold gas in a hot wind" (https://arxiv.org/abs/1806.02728) by Max Gronke and S. Peng Oh.

Absorption and Emission lines of spectra from cold gas around galaxies indicate they possess speeds of about several hunderd km/s. 

Now this causes a trouble from a theoretical point of view which cannot adequately explain the observation. This is because theoretical feedback models like that of Supernova activity can generate hot gas moving at high speed. But the problem is that it cannot efficiently interact with the cold cloud and make it achieve high speeds. This cold gas quickly gets annihilated (crushed) by the hot wind before it can get dragged achieve a high velocity. 

However mixed warm gas (moving with speeds of several hunderd km/s) depending on the prevailing physical condtions can cool efficiently to produce cold gas with high velocity over time which can explain the observation. This is illustrated in the work cited above. Gas clouds of sufficiently large size (few parsecs) can produce even higher amounts of (fast moving) cold gas than what it has started with initially.
