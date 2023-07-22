# Joint caching and routing with any topology: developed poly-time algorithms with approx. guarantees to minimize total cost, formulated as MILP, NP-hard. 

**cite** contains a attack resilience analysis implementation for our [paper](https://nsrg.cse.psu.edu/tbd).  If you find this code useful in your research, please consider citing:

    @INPROCEEDINGS{Xie22ICDCS,
      author = { Tian Xie and Ting He and Patrick McDaniel and Quinn Burke},
      title = {Joint Caching and Routing in Cache Networks with Arbitrary Topology},
      booktitle = {IEEE ICDCS},
      year = {2022},
      month = {July}
    }

This code was tested on Matlab R2019a.


### Structure of Code Design
We have plenty of files in this part, details are explained as following:

we use homogenous cache capacity, which is that all the edge nodes have same storage capacity. All the links initially have the same link capacity. The time unit was set as millisecond.

main three files:

- test_unlimited_link_capacities.m

- test_binary_cache_capacities.m:
    
    - controls link capacity, always keep the two as the same, used to randomize it:
    ```
        linkmin = .2; % .5; % minimum link capacity (before augmentation)
        linkmax = .2; %.5; % maximum link capacity (before augmentation)
        Linkmin = [.2 .4 .6 .8]; %[.5 .75 1 1.25 ];
  ``` 
    - `K_total = [2:2:20]; % optimal K = 12` this is the range of _K_ parameter that you might need to tune. Basically, we should better choose the _K_ value that gives us the minimum congestion. 
    - `K = 20; %12; % optimized for ours` is the selected _K_ that minimize the congestion.
    - `n_cache = 1;` denotes how many nodes are cache.
    - ` cache_pre = 0; % not predefine the cache` is currently not used. Sometimes predetermine which node is the cache. Wasn't sucdessful.
    - `n_client = 5;` is the number of edge nodes excluding the server but also including cache. Should not change.
    
- test_general_case.m: tune the parameter for different network topologies, take `Abovenet` as an example:
    - when you change the `catelog size` eg ` C = 10`, what you should preserve is .
    - these two are control the routing cost as delays in milliseconds. The first one is the local link cost and the second one is the remote link. (make sure that servers are always degree one nodes.)
    ```
        cmin = 1; cmax = 20; % min/max link cost according to [Ioannidis18JSAC], e.g., delay in ms
        cmin_s = 100; cmax_s = 200; % min/max link cost from the remote server (delay in ms)  
  ```
    - this is the skewness, where the reference shows that this is the typical skewness at least for web requests:
    ```
  skewness = 0.7; % typical skewness in web requests [Breslau'99INFOCOM: "Web Caching and Zipf-like Distributions: Evidence and Implications"]  
  ```
    - `c_v = 0;` should be always zero because this is edge caching and all the internal nodes representing backbone POP do not have the cache. All of them are backbone degree zero nodes.
    - `c_client = 2; ` controls non-server edge nodes' cache capacity. The server will always store everything.
    - ` deg_client = 3;` currently not used. Set as upperbound of the degree that served as edge nodes.
    - `n_client = 5;` is the number of edge nodes excluding the server. This number of nodes plus one is the number of nodes with the degree upper bounded by this. (except vario)
    - `k_paths = 10;` is the setting from Ioannidis, the number of candidate path per requester. Either 10 or 30.


### Experiment Procedure
Before starting the experiments, quick mention that you will need to modify the `catalog size`on the edge cache, according to what size you require for your experiment. 

Attention: add matlab_bgl to the path.
```
addpath("matlab_bgl");
```

Our cataloge size right now is 2, too small. at least there should be variations of the three networks. So we will try 10, 20 , 30.
Aim: find the setting that the result is not worse than the previous settings.
The average comparison between ours and the benchmark, and also the absolute congestion level the proposed algorithm can achieve. We don't want to be too much larger than one.

To choose the _K_ parameter: you firstly run the section that gives you the plot of congestion cost vs **_K_** and then find the optimal _K_.

You can run a block of code in matlab where you can type
`Ctrl+Enter` and also use `publish->section break` to choose/create the section you want to run only.
```
%% vary design parameter K:
linkmin = .5; % minimum link capacity (before augmentation)
linkmax = .5; % maximum link capacity (before augmentation)
K_total = [2:2:20]; % optimal K = 12

Cache_pre = setdiff(find(sum(A,2)==1),[server 27 34 36 38 48]);
for preindx = 1:length(Cache_pre)
    cache_pre = Cache_pre(preindx);
    disp(' ')
    disp(['cache = node' num2str(cache_pre) ':'])
```
Run different section for example:
run K tunning experiment-> go back to set the K -> vary the link capacity.

If you change the settings on one, you should change the test setting on other two files.

Duplicate the code folder and do your experiment there in order to always keep the best setting.

Watch out the files' name, every time you run it , data would be overwritten. So please comment that code out or change the filename.

### Plot Explanation
Fig.3 is special case one, unlimited link capacity, proposed should be Alg.1.

Fig.4 porbaly meed to mark the algorithm. Second and third bar are algorithm 2 with different K parameter. The first one is a bound, because that is a LP. In our case is exactly equal to one. The second one is a special case. The third one is a general case, one is varing cache capacity another is varying link capacity. We commented one of them out in our [paper](https://nsrg.cse.psu.edu/files/tbd) since similar result holds. 

Todo 6/24/2021:
work on polish experiments.


### Install Gurobi
Software link: https://www.gurobi.com/downloads/gurobi-software/(or you can use the one I download in the folder '\Gurobi-9.1.2-win64.msi')
Activation license: grbgetkey 22d4083a-de87-11eb-91a2-0242ac130004 
(To install this license on a computer where Gurobi Optimizer is installed, copy and paste the following command to the Start/Run menu (Windows only) or a command/terminal prompt (any system))

Keep in mind the path installed because the remaining set_up work needs to be done in order to start the solver in matlab. 
>Reference: 
>https://www.gurobi.com/documentation/9.1/quickstart_mac/matlab_setting_up_grb_for_.html
>https://www.gurobi.com/documentation/9.1/matlab_html/matlab_setting_up_the_grb_.html

```
>> cd C:\gurobi912\win64\matlab
>> gurobi_setup

The MATLAB interface for Gurobi 9.1.2 has been installed.

The directory
    C:\gurobi912\win64\matlab\
has been added to the MATLAB path.
To use Gurobi regularly, you must save this new path definition.
To do this, type the command
    savepath
at the MATLAB prompt. Please consult the MATLAB documentation
if necessary.
>> cd c:/gurobi912/win64/examples/matlab
>> mip1
          status: 'OPTIMAL'
     versioninfo: [1×1 struct]
         runtime: 0.0020
          objval: 3
               x: [3×1 double]
           slack: [2×1 double]
    poolobjbound: 3
            pool: [1×2 struct]
          mipgap: 0
        objbound: 3
       objboundc: 3
       itercount: 0
    baritercount: 0
       nodecount: 0
x 1
y 0
z 1
Obj: 3.000000e+00

```

In the license, you might need to change the username:
>file:///C:/gurobi912/win64/docs/quickstart/retrieving_a_free_academic.html
```

# DO NOT EDIT THIS FILE except as noted
#
# License ID 662666
# Gurobi license for Pennsylvania State University
ORGANIZATION=Pennsylvania State University
TYPE=ACADEMIC
VERSION=9
HOSTNAME=E5-CSE-364-04
HOSTID=42096ebb
USERNAME=met14_wadmin %tbx5027
EXPIRATION=2022-07-01
KEY=OPBI3X56
CKEY=OPBI3X56
```

# Citation
**cite** contains a attack resilience analysis implementation for our [paper](https://nsrg.cse.psu.edu/files/tbd).  If you find this code useful in your research, please consider citing:

    @INPROCEEDINGS{Xie22ICDCS,
      author = { Tian Xie and Ting He and Patrick McDaniel and Quinn Burke},
      title = {Joint Caching and Routing in Cache Networks with Arbitrary Topology},
      booktitle = {IEEE ICDCS},
      year = {2022},
      month = {July}
    }
