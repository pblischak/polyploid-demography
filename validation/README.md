# Validation Simulations in SLiM

The scripts in this folder implement six main models tested across various
parameter combinations for autotetraploids, allotetraploids, and segmental
allotetraploids. Details on how to run the models can be found in the comments
at the top of each script. In general, the scripts are run by invoking the
``slim`` command from within a terminal and specifying the required parameters
using the ``-d`` flag, followed by the name of the script.

**Example**

```bash
slim -d "rep=1" -d "nuBot=0.1" autotetraploid_bottleneck.slim
```
