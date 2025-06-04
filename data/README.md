# Introduction
For testing the Julia implementation provided within this repository, the following GAP instances from [http://www.co.mi.i.nagoya-u.ac.jp/yagiura/gap/](http://www.co.mi.i.nagoya-u.ac.jp/~yagiura/gap/) are provided: **d05100**, **e201600**, and **e801600**. Several other instances may also be obtained from the aforementioned webpage. Moreover, as described there, data are generating and formatted within each file as follows.

number of agents (m), number of jobs (n)\
for each agent i (i=1,...,m):\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;cost of allocating job j to agent i (j=1,...,n)\
for each agent i (i=1,...,m):\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;resource consumed in allocating job j to agent i (j=1,...,n)\
resource capacity of agent i (i=1,...,m)

Particularly, the first line encodes the number of agents ($i$) and jobs ($j$) found in each data file, while the first loop encodes the appropriate cost ($c_{ij}$) for allocating agent $i$ to job $j$, the second loop encodes the amount of time ($a_{ij}$) used by agent $i$ to complete job $j$, and the final line encodes the available time that agent $i$ has to complete job $j$ within each data file. 

## Pre-processing of GAP Type D Instances 

It is important to state that the type D instance data file contained within this repository (d05100) does not match the file format of those that can be acquired directly from Mutsunori Yagiura's GAP instance website. As such, if one desires to test additional GAP type D instances, then either preprocessing may need to be done on the file so that it matches the format contained in the type E instance files, or a method may need to be written to handle reading in the type D instances. In the latter case, the existing method for reading in data instances (SLBLR/src/get_data.jl) can be modified to handle both type E and type D instances. 
