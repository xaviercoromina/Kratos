import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Settings explanation
# DIMENSION TO COMPUTE:
# This symbolic generator is valid for both 2D and 3D cases. Since the element has been programed with a dimension template in Kratos,
# it is advised to set the dim_to_compute flag as "Both". In this case the generated .cpp file will contain both 2D and 3D implementations.
# LINEARISATION SETTINGS:
# FullNR considers the convective velocity as "v-vmesh", hence v is taken into account in the derivation of the LHS and RHS.
# Picard (a.k.a. QuasiNR) considers the convective velocity as "a", thus it is considered as a constant in the derivation of the LHS and RHS.
# DIVIDE BY RHO:
# If set to true, divides the mass conservation equation by rho in order to have a better conditioned matrix. Otherwise the original form is kept.
# ARTIFICIAL COMPRESSIBILITY:
# If set to true, the time derivative of the density is introduced in the mass conservation equation together with the state equation
# dp/drho=c^2 (being c the sound velocity). Besides, the velocity divergence is not considered to be 0. These assumptions add some extra terms
# to the usual Navier-Stokes equations that act as a weak compressibility controlled by the value of "c".
# CONVECTIVE TERM:
# If set to true, the convective term is taken into account in the calculation of the variational form. This allows generating both
# Navier-Stokes and Stokes elements.

## Symbolic generation settings
mode = "c"
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
ASGS_stabilization = True           # Consider ASGS stabilization terms
formulation = "WeaklyCompressibleNavierStokes" # Element type. Options: "WeaklyCompressibleNavierStokes", "Stokes"

if formulation == "WeaklyCompressibleNavierStokes":
    linearisation = "Picard" # Convective term linearisation type. Options: "Picard", "FullNR"
    output_filename = "weakly_compressible_navier_stokes_volumetric_strain.cpp"
    template_filename = "weakly_compressible_navier_stokes_volumetric_strain_cpp_template.cpp"
else:
    err_msg = f"Wrong formulation. Given '{formulation}'. Available option is 'WeaklyCompressibleNavierStokes'."
    raise Exception(err_msg)

info_msg = "\n"
info_msg += "Element generator settings:\n"
info_msg += "\t - Element type: " + formulation + "\n"
info_msg += "\t - Dimension: " + dim_to_compute + "\n"
info_msg += "\t - ASGS stabilization: " + str(ASGS_stabilization) + "\n"
print(info_msg)

#TODO: DO ALL ELEMENT TYPES FOR N-S TOO
if formulation == "WeaklyCompressibleNavierStokes":
    if (dim_to_compute == "2D"):
        dim_vector = [2]
        nnodes_vector = [3] # tria
    elif (dim_to_compute == "3D"):
        dim_vector = [3]
        nnodes_vector = [4] # tet
    elif (dim_to_compute == "Both"):
        dim_vector = [2, 3]
        nnodes_vector = [3, 4] # tria, tet
elif formulation == "Stokes":
    # all linear elements
    if (dim_to_compute == "2D"):
        dim_vector = [2, 2]
        nnodes_vector = [3, 4] # tria, quad
    elif (dim_to_compute == "3D"):
        dim_vector = [3, 3, 3]
        nnodes_vector = [4, 6, 8] # tet, prism, hex
    elif (dim_to_compute == "Both"):
        dim_vector = [2, 2, 3, 3, 3]
        nnodes_vector = [3, 4, 4, 6, 8] # tria, quad, tet, prism, hex

## Initialize the outstring to be filled with the template .cpp file
print("Reading template file \'"+ template_filename + "\'\n")
templatefile = open(template_filename)
outstring = templatefile.read()

for dim, nnodes in zip(dim_vector, nnodes_vector):

    if(dim == 2):
        strain_size = 3
    elif(dim == 3):
        strain_size = 6

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    ## Unknown fields definition
    v = DefineMatrix('v',nnodes,dim)               # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn',nnodes,dim)             # Previous step velocity
    vnn = DefineMatrix('vnn',nnodes,dim)           # 2 previous step velocity
    eps_vol = DefineVector('eps_vol',nnodes)       # Volumetric strain
    eps_vol_n = DefineVector('eps_vol_n',nnodes)   # Volumetric strain
    eps_vol_nn = DefineVector('eps_vol_nn',nnodes) # Volumetric strain

    ## Fluid properties
    # Note that these are nodal (targeting the two-fluid case)
    k_nodes = DefineVector('k',nnodes)     # Nodal bulk modulus
    c_nodes = DefineVector('c',nnodes)     # Nodal sound speed
    rho_nodes = DefineVector('rho',nnodes) # Nodal density

    k = sympy.Symbol('k', positive = True) # Gauss pt. bulk modulus
    mu = sympy.Symbol('mu', positive = True) # Gauss pt. dynamic viscosity
    rho =sympy.Symbol('rho', positive = True) # Gauss pt. density

    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim)            # Velocity field test function
    q = DefineVector('q',nnodes)                # Pressure field test function

    ## Other data definitions
    f = DefineMatrix('f',nnodes,dim)            # Forcing term

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    h = sympy.Symbol('h', positive = True)             # Element size
    dt  = sympy.Symbol('dt', positive = True)          # Time increment
    dyn_tau = sympy.Symbol('dyn_tau', positive = True) # Stabilization dynamic tau constant
    stab_c1 = sympy.Symbol('stab_c1', positive = True) # Stabilization constant 1
    stab_c2 = sympy.Symbol('stab_c2', positive = True) # Stabilization constant 2

    ## BDF time discretization coefficients
    bdf0 = sympy.Symbol('bdf0')
    bdf1 = sympy.Symbol('bdf1')
    bdf2 = sympy.Symbol('bdf2')

    ## Data interpolation to the Gauss points
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N
    eps_vol_gauss = eps_vol.transpose()*N

    ## Convective velocity definition
    if linearisation == "Picard":
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol
    elif (linearisation == "FullNR"):
        vmesh = DefineMatrix('vmesh',nnodes,dim)    # Mesh velocity
        vconv = v - vmesh                           # Convective velocity defined as a velocity dependent variable
    else:
        raise Exception("Wrong linearisation \'" + linearisation + "\' selected. Available options are \'Picard\' and \'FullNR\'.")
    vconv_gauss = vconv.transpose()*N

    ## Compute the stabilization parameters
    stab_norm_a = 0.0
    for i in range(0, dim):
        stab_norm_a += vconv_gauss[i]**2
    stab_norm_a = sympy.sqrt(stab_norm_a)
    tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*stab_norm_a)/h + (stab_c1*mu)/(h*h)) # Stabilization parameter 1
    # tau2 = (mu + (stab_c2*rho*stab_norm_a*h)/stab_c1)                                # Stabilization parameter 2
    tau2 = 0.0 #FIXME: Design a tau for the volumetric strain

    ## Compute the accelerations at Gauss points
    accel_gauss = (bdf0*v + bdf1*vn + bdf2*vnn).transpose()*N
    accel_eps_vol_gauss = (bdf0*eps_vol + bdf1*eps_vol_n + bdf2*eps_vol_nn).transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    grad_w = DfjDxi(DN,w)
    grad_q = DfjDxi(DN,q)
    grad_v = DfjDxi(DN,v)
    grad_eps_vol = DfjDxi(DN,eps_vol)

    div_w = div(DN,w)
    div_v = div(DN,v)
    div_vconv = div(DN,vconv)

    grad_sym_v = grad_sym_voigtform(DN,v)       # Symmetric gradient of v in Voigt notation
    grad_w_voigt = grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation
    # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.

    # Convective term definition
    convective_term_gauss = vconv_gauss.transpose()*grad_v

    ## Compute galerkin functional
    # Navier-Stokes functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss
    rv_galerkin -= rho*w_gauss.transpose()*accel_gauss
    rv_galerkin -= rho*w_gauss.transpose()*convective_term_gauss.transpose()
    rv_galerkin -= grad_w_voigt.transpose()*stress
    rv_galerkin -= div_w*k*eps_vol_gauss
    rv_galerkin -= q_gauss*accel_eps_vol_gauss
    rv_galerkin += q_gauss*div_v

    ##  Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used
    mom_residual = rho*f_gauss
    mom_residual -= rho*accel_gauss
    mom_residual -= rho*convective_term_gauss.transpose()
    mom_residual += k*grad_eps_vol

    # Mass conservation residual
    mass_residual = accel_eps_vol_gauss
    mass_residual -= div_v

    vel_subscale = tau1*mom_residual
    eps_subscale = tau2*mass_residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    rv_stab = rho*div_vconv*w_gauss.transpose()*vel_subscale
    rv_stab += rho*vconv_gauss.transpose()*grad_w*vel_subscale
    rv_stab -= div_w*k*eps_subscale
    rv_stab -= grad_q.transpose()*vel_subscale

    ## Add the stabilization terms to the original residual terms
    if ASGS_stabilization:
        rv = rv_galerkin + rv_stab
    else:
        rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes*(dim+1), 1)
    testfunc = sympy.zeros(nnodes*(dim+1), 1)

    for i in range(0,nnodes):
        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]
        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = eps_vol[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    print("Computing " + str(dim) + "D RHS Gauss point contribution\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    print("Computing " + str(dim) + "D LHS Gauss point contribution\n")
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace("//substitute_lhs_" + str(dim) + 'D' + str(nnodes) + 'N', lhs_out)
    outstring = outstring.replace("//substitute_rhs_" + str(dim) + 'D' + str(nnodes) + 'N', rhs_out)

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()
