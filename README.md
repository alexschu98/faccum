# faccum

faccum is a tiny library I wrote during my master's thesis to calculate the generalized factorial cumulants of a stochastic system.

## Installation
Clone the repo and import the file.

## Usage 

Create your W_z matrix dependent on z as a function and just call the function.

```python
import faccum

# setup w matrix
def Wz(z): 
    ...

# get faccums up to m = 5 at t = 100
fc = faccum.generalized_factorial_cumulants(Wz, t=100, s0=1, mmax=5)
```

 ## License
[MIT License](https://opensource.org/licenses/MIT)