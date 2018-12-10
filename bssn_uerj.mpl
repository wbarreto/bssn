##############################################################
# BSSN equations in spherical symmetriy for the EKG system,  #
# by W. Barreto, UERJ 16/10/2018.                            #
# Based on the Maple script of Akbarian & Choptuik           #
# for the paper PRD 121, 021103 (2015),                      #
# whom calculate using GRTENSOR.                             #
##############################################################

with(tensor):

# Definitions 

coord := [x, theta, phi]:

alias(A = A(t, x)):        # 1 conformal 3-metric ; \tilde \gamma_{rr}
alias(B = B(t, x)):        # 2 conformal 3-metric
alias(Lam = Lam(t,x)):     # 3 conformal connection 
alias(Gam = Gam(t,x)):     # connection; this is not used indeed! Only for testing
alias(Phi = Phi(t,x)):     # 4 Conformal factor
alias(Ka = Ka(t,x)):         # 5 Extrinsic Curvature
alias(Axx = Axx(t,x)):     # 6 Conformal Extrinsic Curvature 
alias(Athth = Athth(t,x)): # 7 Conformal Extrinsic Curvature 

alias(beta = beta(t,x)):   # Shift
alias(alpha = alpha(t,x)): # Lapse

alias(S = S(t,x)):         # Source 

alias(divbeta = divbeta(t,x)): # Only as temporal variable

####################
# Conformal metric #
####################

g_compts := array(symmetric, sparse, 1 .. 3, 1 .. 3):

g_compts[1,1] := A:
g_compts[2,2] := x^2*B:
g_compts[3,3] := x^2*sin(theta)^2*B:

g := create([-1, -1], eval(g_compts)):

ginv := invert(g, 'detg'):

D1g := d1metric ( g, coord ):

D2g := d2metric ( D1g, coord ):

Cf1 := Christoffel1 ( D1g ):

Cf2 := Christoffel2 (ginv,Cf1):

################
# Flat metric  #
################
gf_compts := array(symmetric, sparse, 1 .. 3, 1 .. 3):

gf_compts[1,1] := 1:
gf_compts[2,2] := x^2:
gf_compts[3,3] := x^2*sin(theta)^2:

gf := create([-1, -1], eval(gf_compts)):

gfinv := invert(gf, 'detgf'):

D1gf := d1metric ( gf, coord ):

D2gf := d2metric ( D1gf, coord ):

Cf1f := Christoffel1 ( D1gf ):

Cf2f := Christoffel2 (gfinv,Cf1f):

#----------------------------------

# Definitions of \tilde\Gamma and \tilde\Lambda
# [Eqs. (11) and (12) in AC paper ]

gCf2:=prod(ginv,Cf2):

TG:=contract(gCf2,[1,4],[2,5]):

gCf2f:=prod(ginv,Cf2f):
TGf:=contract(gCf2f,[1,4],[2,5]):

TG2:=lin_com(TG,-1,TGf):

Lamb_compts := array(sparse,1..3):
Lamb_compts[1]:=Lam:
Lamb := create([1], eval(Lamb_compts)):

Gamm:=lin_com(Lamb,TGf):

###################################################################
# Calculation of the 3-Ricci tensor decomposed                    #
# as R_{ij}=R^\Phi_{ij} + \tilde R_{ij} [Eq. (21), (22) and (23)] #
###################################################################

# Here be careful, because indexing of Cf1 is not the
# same of Akbarian & Choptuik [Cf1_{ijk} -> ChL_{kij}]

CC1:=prod(Cf2,Cf1):
CC2:=prod(ginv,CC1):
CC3:=symmetrize(CC2,[5,8]):
CC4:=contract(CC3,[1,4],[2,7],[3,6]):
CC5:=contract(CC2,[1,6],[2,5],[3,8]):

SS:=lin_com(2,CC4,1,CC5):
#ss:=get_compts(SS):

CC6:=symmetrize(Cf1,[1,3]):
CC7:=prod(Gamm,CC6): 
CC7:=contract(CC7,[1,3]):

CC8:=partial_diff(Gamm,coord):
CC9:=prod(g,CC8):
CC10:=symmetrize(CC9,[2,4]):
CC11:=contract(CC10,[1,3]):

CC12:=partial_diff(g,coord):
CC13:=partial_diff(CC12,coord):
CC14:=prod(ginv,CC13):
CC15:=contract(CC14,[1,5],[2,6]):

# [Eq. (23)]

RB:=lin_com(-1/2,CC15,1,CC11,1,CC7,1,SS):

rb:=get_compts(RB):

#expand(rb[1,1]);
#expand(rb[2,2]);

CF:=create([],Phi):
DCF:=cov_diff(CF,coord,Cf2):

CC16:=cov_diff(DCF,coord,Cf2):
CC17:=raise(ginv,CC16,1):
CC18:=contract(CC17,[1,2]):
CC19:=prod(g,CC18):
CC20:=raise(ginv,DCF,1):
CC21:=prod(CC20,DCF):
CC22:=contract(CC21,[1,2]):
CC23:=prod(DCF,DCF):
CC24:=prod(g,CC22):

# [Eq.(22)]

RPhi:=lin_com(-2,CC16,-2,CC19,4,CC23,-4,CC24):

rphi:=get_compts(RPhi):

#expand(rphi[1,1]);
#expand(rphi[2,2]);
#expand(rphi[3,3]);

# [Eq. (23)] 

# Here AC calculate R_{ij} to make a test...

########################### 
# Divergence of the Shift #
###########################

SH_compts := array(sparse,1..3):
SH_compts[1] := beta:
SH := create([1],eval(SH_compts)):

DSH:=cov_diff(SH,coord,Cf2):
DIVSH:=contract(DSH,[1,2]):

DIVBETA:=create([],divbeta):

# Matter field [Eq. (25)]

S_compts := array(sparse, 1..3):
S_compts[1] := Sx:
S := create([1],eval(S_compts)):

# Traceless conformal (covariant) extrinsic curvature

A_compts:=array(symmetric,sparse,1..3,1..3):

A_compts[1,1]:=Axx:
A_compts[2,2]:=x^2*Athth:
A_compts[3,3]:=x^2*Athth*sin(\theta)^2:

A := create([-1,-1],eval(A_compts)):

Aud := raise(ginv,A,1);
Auu := raise(ginv,Aud,2);

# Lapse

Lapse:=create([],alpha):

# K

K:=create([],Ka):

#####################################################
# First BSSN equation [Eq. (18) in the paper of AC] #
# (BSSN Eq. for \tilde\Gamma) --> EE1 (below)       #
#####################################################

EE11:=partial_diff(SH,coord):
EE12:=partial_diff(EE11,coord):
EE13:=prod(ginv,EE12):

TEE1:=contract(EE13,[1,4],[2,5]):

EE21:=partial_diff(Lapse,coord);
EE22:=prod(Auu,EE21):

TEE2:=contract(EE22,[2,3]):
Cf2;
EE31:=prod(Cf2,Auu);
EE32:=contract(EE31,[3,4],[2,5]):

EE33:=partial_diff(K,coord):
EE34:=prod(ginv,EE33):
EE35:=conract(EE34,[2,3]):

EE36:=prod(ginv,S):
EE37:=contract(EE36,[2,3]):

EE38:=prod(A,DCF):
EE39:=contract(EE38,[2,3]):

EE310:=lin_com(1,EE32,-2/3,EE35,-1,EE37,6,EE39):

TEE3:=prod(Lapse,EE310):

EE41:= prod(SH,Gamm):

TEE4:= contract(EE41,[1,3]):

EE51:=prod(Gamm,SH):

TEE5:=contract(EE51,[1,3]):

EE61:=prod(Gamm,DIVBETA):
EE62:=partial_diff(DIVBETA,coord):
EE63:=prod(ginv,EE62):
EE64:=contract(EE63,[2,3]):

TEE6:=lin_com(2,EE61,1,EE64):

# MCh{^i j k}*A{^k ^j} - 2/3*g{^i ^j}*Extc{,j} - g{^i ^j}*S{j} + 6*A{^i ^j}*DCF{j}

# Sh{^j}*GA{^i ,j}

# GA{^j}*Sh{^i ,j}

# 2*GA{^i}*divbeta(t,x)+g{^i ^l}*DBS{,l})

EE1:=lin_com(1,TEE1,-2,TEE2,2,TEE3,1,TEE4,-1,TEE5,1/3,TEE6):

