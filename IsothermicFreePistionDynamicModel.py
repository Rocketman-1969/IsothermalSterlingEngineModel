# Isothermic free pistion dynamic model
#alsdfjasl;fdj
#imports 
#Reads json input files
import json
#basic math and matrix operations
import numpy as np
#Complex math ops ODE solver and integrator
import scipy
#plotting tools
from matplotlib import pyplot as plt

#initilialize class
class DyFPSolv:
    
	#pulls Data from input file
    def __init__(self,filename):
         
         #read json input file
        json_data = json.loads(open(filename).read())

#---------Temperature Data [K]----------------------------------------

        #Hot side
        self.T_h=json_data["T_h[k]"]

        #Cold Side
        self.T_k=json_data["T_k[k]"]

        #Regenerator
        self.T_r=(self.T_h-self.T_k)/(np.log(self.T_h/self.T_k))
        print(self.T_r)
#---------Pressure Data [Bar]----------------------------------------

        #Average Pressure
        self.P_avg=json_data["P_avg[Bar]"]

#---------Regenerator Porosity [%]-----------------------------------

        #Porosity of Regenerator Mesh
        self.Porosity=json_data["Porosity[%]"]

#---------Engine Diameter Data [cm]----------------------------------

        #Diameter of Power Pistion
        self.d_p=json_data["d_p[cm]"]

        #Diameter of Displacer Pistion
        self.d_d=json_data["d_d[cm]"]

        #Diameter of Displacer Rod
        self.d_rod=json_data["d_r[cm]"]

        #Diameter of Hot Tube
        self.d_h=json_data["d_h[cm]"]

        #Diameter of Regenerator Wire Mesh
        self.d_w=json_data["d_w[cm]"]

#---------Engine Length Data [cm]-----------------------------------

        #Length of Cold Tube
        self.L_k=json_data["L_k[cm]"]

        #Length of Hot Tube
        self.L_h=json_data["L_h[cm]"]

        #Length of Regenerator
        self.L_r=json_data["L_r[cm]"]

        #Clearance Distance in Expansion Space
        self.C_e=json_data["C_e[cm]"]

        #Clearance Distance in Compression Space
        self.C_c=json_data["C_c[cm]"]

        #Maximum Displacement of Power Piston
        self.X_p=json_data["X_p[cm]"]

        #Maximum Displacement of Displacer 
        self.X_d=json_data["X_p[cm]"]

#---------Engine Area Data [cm^2]-----------------------------------

        #Cross Sectional Area of Hot Tube
        self.A_h=json_data["A_h[cm2]"]

        #Cross Sectional Area of Regenerator
        self.A_r=json_data["A_r[cm2]"]

        #Cross Sectional Area of Power Pistion
        self.A_p=np.pi/4*self.d_p**2

        #Cross Sectional Area of Displacer
        self.A_d=np.pi/4*self.d_d**2

        #Cross Sectional Area of Rod
        self.A_rod=np.pi/4*self.d_rod**2

        #Cross Sectiional Area of cold tube
        self.A_k=np.pi/4*self.d_h**2

        #Wetted Area of Cold Tube
        self.A_W_k=json_data["A_W_k[cm2]"]

        #wetted area of hot tube
        self.A_W_h=np.pi*self.d_h*self.L_h

#---------Engine Volume Data [cm^3]---------------------------------

        #Average Volume of Bounce Space
        self.V_B=json_data["V_B[cm3]"]

        #Average Volume of Dead Space (gas spring)
        self.V_D=json_data["V_D[cm3]"]

        #Volume of Regenerator
        self.V_r=json_data["V_r[cm3]"]

        #Volume of the cold tube
        self.V_k=self.A_h*self.L_k

        #Volume of the hot tube
        self.V_h=self.A_h*self.L_h


#---------Operation Perameters--------------------------------------

        #Load Damping[Ns/m]
        self.F_L=json_data["F_L[Ns/m]"]

        #Phase Angle between Power Piston and Displacer[Degrees]
        self.phi=json_data["Phi[degrees]"]

        #Engine Operating Freq [rad/s]
        self.omega=json_data["omega[rad/s]"]

        #Engine Working Gas
        self.WorkingGas=json_data["WorkingGas"]

        #Dynamic Viscocity at T_k[kg/m/s]
        self.nuT_k=json_data["nuT_k"]

        #Dynamic Viscocity at T_h[kg/m/s]
        self.nuT_h=json_data["nuT_h"]

        #Dynamic Viscocity at T_r[kg/m/s]
        self.nuT_r=json_data["nuT_r"]

        #Density of helium at 50 bar T_k[kg/m3]
        self.rhoT_k=json_data["RhoT_k"]

        #density of helium at 50 bar T_h[kg/m3]
        self.rhoT_h=json_data["RhoT_h"]

        #Density of helium at 50 bar T_r[kg/m3]
        self.rhoT_r=json_data["RhoT_r"]



#--------Wetted Diameters------------------------------------------
        
        #wetted diameter of cold tube
        self.D_H_k=(self.V_k*4)/self.A_W_k

        #wetted diameter of hot tube
        self.D_H_h=(self.V_h*4)/self.A_W_h

        #wetted diameter of regenerator
        self.D_H_r=(self.d_w*self.Porosity*.01)/(1-self.Porosity*.01)


    #Calculates maximum Velocity of the gas throught hot, cold, and regenerator
    def MaxVelocity(self):

        #equation is broken down for readablity
        #parts are arbitrarly assigned A, B, and C
        A=(self.A_p*self.X_p)**2
        B=((2*self.A_d-self.A_rod)*self.X_d)**2
        C=2*self.A_p*self.X_p*(2*self.A_d-self.A_rod)\
            *self.X_d*np.sin(np.deg2rad(self.phi))
         
        #Maximum velocity through hot tube
        self.U_h=(self.omega*np.sqrt(A+B-C))/self.A_h

        #Maximium velocity through regenerator
        self.U_r=(self.omega*np.sqrt(A+B-C))/self.A_r

        #Maximum velocity through cold tube
        self.U_k=(self.omega*np.sqrt(A+B-C))/self.A_k

    #Calculates the instantanious volume flow rate from expanstion to compression
    def InstVolFloRate(self, xdot_p, xdot_d):

        Vdot=self.A_p*xdot_p-(2*self.A_d-self.A_rod)*xdot_d

        return Vdot
    
    #calculates the instantanious velcoity of the gas through hot, cold and regenerator
    def InstVelocity(self,Vdot):

        #inst velocity through hot tube
        u_h=Vdot/self.A_h

        #inst velocity through cold tube
        u_k=Vdot/self.A_k

        #inst velocity through regenerator
        u_r=Vdot/self.A_r

        return u_h, u_k, u_r

    #calculates Darcy friction factor
    def C_f(self,u_h,u_k,u_r):
        
        #calculate Reynolds number for hot tube 
        #dynamic viscosity at T_h
        Re_h=(u_h*self.L_h*.01)/self.nuT_h

        #calculate Reynolds number for cold tube
        #dynamic viscosity at T_k
        Re_k=(u_k*self.L_k*.01)/self.nuT_k

        #calculate reynolds number for regenerator
        #dynamic viscosity at T_r
        Re_r=(u_r*self.L_r)/self.nuT_r

        #Calculate Darcy Friction factor based on Reynolds number
        #hot tube
        if Re_h<=2000: 
            Cf_h=64/Re_h
        else:
            Cf_h=0.316*Re_h**(-.25)

        #cold tube
        if Re_k<2000: 
            Cf_k=64/Re_k
        else:
            Cf_k=0.316*Re_k**(-.25)

        #Regeneratior
        if Re_r<=60:
            Cf_r=4*10**(1.73-.93*np.log10(Re_r))

        elif 60<Re_r<=1000:
            Cf_r=4*10**(.714-.365*np.log10(Re_r))

        else:
            Cf_r=4*10**(.0015-.125*np.log10(Re_r))
        
        return Cf_h, Cf_k, Cf_r

    #calculate the total pressure drop in Pa
    def DeltaP(self,Cf_h,Cf_k,Cf_r,u_k,u_r,u,h):
        P_k=.01*.5*self.rhoT_k*((Cf_k*self.L_k)/self.D_H_k)*np.abs(u_k)*u_k

        
        
        return 0
        

    def test(self):
        return self.T_r


if __name__=="__main__":
    solver=DyFPSolv("inputfile.json")
    T_r=DyFPSolv.test
    print(T_r)
	