
file="cantilever.msh";

# MATERIAL: Young and Poisson
material = [[1000, 0., 0.1]];

# SOLID: associates materials to sets
solid = [[4, 1]]

# nodal DBC: each row a bc. node, direction, val
dbc = [[1,1],[1,2]]  
dbcval = [0.,0.]
S0=[0.0; 0.0; 0.0]  

# eleset where follower forces are applied
tbc = [3]
# define parameter for the analysis
info=Sinfo()

# Mass damped - 1%
# Pcoal = 7.12118403
# Pbif = 7.122053
# info.alpha = 0.68/50   # Mass damping coefficient
# info.beta = 0          # Stiffness damping coefficient

# Mass damped high - 20%
Pcoal = 7.12118403
Pbif = 7.45455383
info.alpha = 0.68/2.5  # Mass damping coefficient
info.beta = 0          # Stiffness damping coefficient

# Stiffness damped - 10%
# Pcoal = 7.21373583
# Pbif = 7.21373583
# info.alpha = 0           # Mass damping coefficient
# info.beta = 1.0/(0.68*5) # Stiffness damping coefficient

# tbcval = [0.95*Pcoal] # Before
# tbcval = [Pcoal]      # CP
tbcval = [Pbif]         # BP
# tbcval = [1.05*Pbif]  # After
# tbcval = [1.15*Pbif]  # After far
# tbcval = [1.3*Pbif]   # 1.3 P
tbcval_damping = [Pbif] # Load value for which damping is computed. Only pertinent for stiffness damping

info.neig=2    # number of modes to be computed
info.Lmm = [1,2]  # list of the computed modes retained as master
info.Ffreq=1  # mode number that will give the freq (only one but +iomegat and -iomegat)
#info.Fmodes=[1]  # loading will be prop to sum of these modes  
#info.Fmult=[0.5]  # with these amplitudes
info.style = 'c'
info.max_order = 9
info.max_orderNA = 1
info.tol=1e-1

info.LambdaSym = [1,1,-1,-1,0]
# info.LambdaSym = [1,-1,0] # when info.Lmm = [1]
info.Jordan = true
info.symbolic_Jordan = true                       # Defines wheather a numerical or symbolic check of Jordan blocks is performed
info.tolHopf = 10^-2
info.Hopf_pairs = [1 3; 2 4]                      # Only needed if Jordan check is symbolic