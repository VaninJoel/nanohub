
from cc3d.core.PySteppables import *
from numpy.random import randint

# Antimony String Model reference for the Tellurium notebook 
# (TODO: add shared *url*)
model_pyro = '''
    # Equations     
    J2_3_1:  -> NLRP3i;  a1 * (((NF_kBn - NF_kB0)^gamma_NF) / (NF_50^gamma_NF + ((NF_kBn - NF_kB0)^gamma_NF)));
    J2_3_2:  NLRP3i -> ;       delta_1 * NLRP3i;
    J2_3_3:  NLRP3i -> NLRP3a; S2 * k1 * NLRP3i;
    
    J2_4_1:  NLRP3a -> ;       delta_1 * NLRP3a;
    J2_4_2:  NLRP3a -> NLRP3o; k2 * NLRP3a^2;  
    

    J2_8:    ASCf -> ASCb;     k3 * (1 / (1+((NLRP3o + a)/ b )^-c)) * NLRP3o * ASCf; 

   
    J2_13:   pro_C1 -> C1;     k4 * ASCb * pro_C1;      


    J2_16:   GSDMD -> GSDMD_N; a2 * ((C1^gamma_C1) / ((C1_50^gamma_C1) + (C1^gamma_C1))) * GSDMD;
    
    
    J2_18_1: -> pro_IL_1beta;  a3 * ((NF_kBn - NF_kB0)^gamma_NF / (NF_50^gamma_NF + (NF_kBn - NF_kB0)^gamma_NF));
    J2_18_2: pro_IL_1beta -> ; delta_2 * pro_IL_1beta;        
    J2_18_3: pro_IL_1beta -> IL_1betac; a4 * ((C1^gamma_C1) / ((C1_50^gamma_C1) + (C1^gamma_C1))) * pro_IL_1beta;    
    
    
    J2_19:   IL_1betac -> ;             delta_2 * IL_1betac;    
    J2_20:   IL_1betac -> IL_1betae;    k5 * (GSDMD_N / (GSDMD + GSDMD_N)) * IL_1betac;
        
    
    J2_23:   pro_IL_18 -> IL_18c;       a5 * ((C1^gamma_C1) / ((C1_50^gamma_C1) + (C1^gamma_C1))) * pro_IL_18;    
    
    J2_24:   IL_18c -> IL_18e;          k6 * (GSDMD_N / (GSDMD + GSDMD_N)) * IL_18c;
    
    
    J2_26:   -> Volume;                 k7 * (GSDMD_N / (GSDMD + GSDMD_N)) * Volume;

    # Constants
    a = 1;           # arbitrary units (a.u.)
    b = 2;           # (a.u.)
    c = 1000;      
    
    a1 = 0.07;       # min^-1 (a.u.)
    a2 = 0.1;        # min^-1
    a3 = 0.06;       # min^-1 (a.u.)
    a4 = 1;          # min^-1
    a5 = 1;          # min^-1

    C1_50 = 0.3;     # (a.u.)
    NF_50 = 0.3;     # (a.u.)

    delta_1 = 0.002; # min^-1
    delta_2 = 0.004; # min^-1

    h = 0.55;        # (a.u.)

    gamma_C1 = 2;    
    gamma_NF = 2; 

    k1 = 0.7;        # min^-1
    k2 = 1;          # min^-1 (a.u.)
    k3 = 0.04;       # min^-1 (a.u.)
    k4 = 0.03;       # min^-1 (a.u.)
    k5 = 1;          # min^-1
    k6 = 1;          # min^-1
    k7 = 0.2;        # min^-1

    n = 1;           # (a.u.)    

    s = 0.8;        # (a.u.) paper
    
    tau = 10;        # min
    Vc = 1.5;        # (a.u.)    

    # Equations 
    NF_kBn  := NF_kB0 +  S1 * h * expo;
    
    log_two := log((T_time / tau));
    expo    := exp(-(log_two^2 / s));
    T_time  := time+1;
        
    # Initial Parameters
    S1        = 1;
    S2        = 1;
    NF_kB0    = 0.25;   # (a.u.) 
    GSDMD     = 1;      # (a.u.)
    ASCf      = 1;      # (a.u.)
    pro_C1    = 1;      # (a.u.)
    Volume    = 1;      # (a.u.)    
    pro_IL_18 = 1;      # (a.u.)
 
    
'''

# Threshold for Volume, ILbeta, DAMP defined as CONSTANTS 
CRITICAL_VOLUME = 1.5 * (25)    # (A.U. from paper) (25 + (25/2)) (1.5 = +50%)  
TLRI_CRITICAL = 0.1             # (A.U.) TODO: NOT CALIBRATED should TIME matter and not only CONCENTRATION to calculate this value?
TLRD_CRITICAL = 0.1             # (A.U.) TODO: NOT CALIBRATED

VIRAL_BURST = 1.0               # (A.U.) Can be changed by user using slider
VIRAL_PRODUCTION = 0.005        # (A.U.) Can be changed by user using slider

class celldeath_v3Steppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
    def diffusionDetection(self, cell, IL1beta, DAMP): 
        """
        Detects amount of IL1beta and DAMP seen by cell in its center each MCS.
        
        Parameters
        ----------
        cell : self.cell object
            A CC3D cell object
        IL1beta : self.field object
            A CC3D field of diffusion object
        DAMP : self.field object
            A CC3D field of diffusion object

        Returns
        -------
        list
            with the values of [IL1beta, DAMP] seen by cell
        """                 
        return [IL1beta[int(cell.xCOM), int(cell.yCOM), 0],
                DAMP[int(cell.xCOM), int(cell.yCOM), 0]] 

    def start(self):
        """
        initialization of diffusion fields, cells and plots     
               
        Yields
        ------
        self.field object
            DAMP is representing viral load in this case the 
            name was inherited from (Hamis & Macfarlane, 2021)
        self.cell object
            all cells start as healthy cells (type=1) from XML
        self.plot_win object
            tracks the number of cells of each type at mcs
        """
                 
        # Diffusion fields       
        DAMP = self.field.DAMP # "virus"
        
        # Initiating Diffusion position and amount        
        DAMP[85, 50, 0] = 30   # center-right
        # DAMP[30, 85, 0] = 30 # upper-left
        # DAMP[30, 50, 0] = 30 # lower-centerish       
        

        for cell in self.cell_list:         
           
           cell.targetVolume = 25  # (A.U.) for spacial
           cell.lambdaVolume = 2.0
        
        # Plot initiation
        self.plot_win = self.add_new_plot_window(title='Cell count',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True,
                                                 config_options={'legend': True})

        self.plot_win.add_plot("Healthy", style='Lines', color='blue', size=1)
        self.plot_win.add_plot("Pyro_DAMP", style='Lines', color='green', size=1)
        self.plot_win.add_plot("Pyro_IL1beta", style='Lines', color='red', size=1)
        self.plot_win.add_plot("Dead_Pyro", style='Lines', color='yellow', size=1)

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        
        # Antimony Model Step
        self.timestep_sbml()
        
        # Diffusion Values
        IL1beta = self.field.IL1B        
        DAMP = self.field.DAMP
        
        
        
        
        
        # Secretor Definition
        secretor_DAMP = self.get_field_secretor("DAMP")
        secretor_IL1beta = self.get_field_secretor("IL1B")
        
        # Plot data
        cell_type_plot = [0, 0, 0, 0]        
        
        for cell in self.cell_list:            
            
            # Possible TODO Random event to initiate S2 instead of assuming S2 = 1
            # Random number generator (1, 100)
            # cell.dict['RandEvent'] = randint(1, 100)                     
            
            # ---- Healthy cell steps ----
            if cell.type == 1: # "Healthy" 

                # Plot information
                cell_type_plot[0] += 1
      
                # Check diffusion field presence in the cell
                il1b_damp = self.diffusionDetection(cell, IL1beta, DAMP)
                
                cell.dict['TLRi'] = il1b_damp[0] # TLR triggered by IL1beta (IL1B)
                   
                cell.dict['TLRd'] = il1b_damp[1] # TLR triggered by DAMPS ("Virus")
                
                # S2 binary function TODO: Implement the new S2 interaction
                # if cell.dict['RandEvent'] > 50: # TODO define percentage
                    # cell.dict['S2'] = 1
                
                # Check for cell type change
                if cell.dict['TLRi'] > TLRI_CRITICAL:  # TODO define resonable (A.U.)          
                    cell.type = 3 # Pyro_IL1Beta
                    
                    # Loading the Pyro_model into the cell
                    self.add_antimony_to_cell(model_string=model_pyro, model_name='m',
                                                cell=cell, step_size=1.0)
                    
                    cell.targetVolume = cell.sbml.m['Volume'] * 25  # 25(A.U.) for spacial representation
                    cell.lambdaVolume = 2.0 # (A.U.) for spacial representation 
                    
                    # Depricated since antimony: cell = PyroptosisCell.cell_init(self,cell)          
                
                elif cell.dict['TLRd'] > TLRD_CRITICAL: # TODO define resonable (A.U.)
                    cell.type = 2 # Pyro_DAMPS 
                    
                    # Loading the Pyro_model into the cell    
                    self.add_antimony_to_cell(model_string=model_pyro, model_name='m',
                                                cell=cell, step_size=1.0)
                                                
                    cell.targetVolume = cell.sbml.m['Volume'] * 25  # 25(A.U.) for spacial representation
                    cell.lambdaVolume = 2.0 # (A.U.) for spacial representation 
          
          
            # ---- Pyro_IL1beta and Pyro_DAMPS cell steps ----         
            elif cell.type == 2 or cell.type == 3:
                
                # Plot information
                if cell.type == 2:
                    cell_type_plot[1] += 1
                else:
                    cell_type_plot[2] += 1                    
                
                # Check for cell type change
                if cell.sbml.m['Volume'] > CRITICAL_VOLUME: # Volume is responsible for death 
                    
                    if cell.type == 2:                        
                        # Simulate a burst of "Virus" DAMPS on cells that were initiated Pyro by DAMPS
                        secretor_DAMP.secreteOutsideCellAtBoundary(cell, VIRAL_BURST)
                    
                    cell.type = 4 # "DEAD_PYRO"
                
                else:
                    cell.targetVolume = 25 * cell.sbml.m['Volume']
                
                # Cell Diffusion Update
                secretor_DAMP.secreteOutsideCellAtBoundary(cell, VIRAL_PRODUCTION)
                secretor_IL1beta.secreteOutsideCellAtBoundary(cell, cell.sbml.m['IL_1betae'])
                
                
            # ---- Dead cells steps ---- 
            else:
                # TODO need to develop the result of a DEAD_PYRO
                # What kind of inflamossomes are released in the external enviroment
                # NEED another type of diffusion maybe?

                # Plot information
                cell_type_plot[3] += 1
               
                
        # ---- Plot Step ----       
        self.plot_win.add_data_point("Healthy",      mcs, cell_type_plot[0])
        self.plot_win.add_data_point("Pyro_DAMP",    mcs, cell_type_plot[1])
        self.plot_win.add_data_point("Pyro_IL1beta", mcs, cell_type_plot[2])
        self.plot_win.add_data_point("Dead_Pyro",    mcs, cell_type_plot[3])      
    
    
    def add_steering_panel(self):
        
        self.add_steering_param(name='Viral Release at Cell Death', val=VIRAL_BURST, min_val=0.1, max_val=50.0,
                                decimal_precision=2, widget_name='slider')
        self.add_steering_param(name='Viral Production at Infected Cell', val=VIRAL_PRODUCTION, min_val=0.001, max_val=0.1,
                                decimal_precision=3, widget_name='slider')
        self.add_steering_param(name='IL1B Decay', val=0.0, min_val=0.0, max_val=0.01,
                                decimal_precision=4, widget_name='slider')
    
    def process_steering_panel_data(self): 
        
        global VIRAL_PRODUCTION, VIRAL_BURST;
        
        VIRAL_PRODUCTION = self.get_steering_param('Viral Production at Infected Cell')
        VIRAL_BURST =      self.get_steering_param('Viral Release at Cell Death')
        decay = self.get_xml_element('il1b_decay')
        decay.cdata = self.get_steering_param('IL1B Decay')
   
    def finish(self):
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return



