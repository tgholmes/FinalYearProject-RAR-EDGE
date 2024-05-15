# Radial-Acceleration-Relation of EDGE simulations

# This project was my final year project (dissertation) for my undergraduate physics degree at the University of Surrey (2024), supervised by Justin Read
# EDGE simulations are from cross-university project, consisting of researchers Justin Read, Oscar Agertz and more
# High-resolution simulations of dwarf galaxys have been produced recorded (see Matthew Orkney et al, MNRAS, 504 (2020))
# The Radial Acceleration Relation describes relationship between the baryonic (expected) gravitational acceleration and the observed gravitational acceleration within galaxies
# Larger (real) galaxies have been found to follow a strong relationship, and we wish to see if our dwarf simulations align with this relation

# Within this main code file, we:
# - Read in the .data files for the galaxy simulations, assort the data into relevant dictionaries and classes
# - Visualise the data using X,Y plots of stars within the simulations
# - Calculate the baryonic and observed accelerations using the radii-dependent cumulative mass for the galaxies

# The EDGE simulations have a LoRes and HiRes counterpart. The LoRes files were small enough to be calculated and compiled within this code
# The HiRes files were too large and too computationally intensive, so they underwent a data refining/reduction process that can be seen elsewhere on the project GITHUB
# These files had their accelerations calculated seperately

#---------------------------------------------------------------------------------

# Importing relevant modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Editing font for all plots produced (TNR chosen due to use within scientific report)
plt.rcParams["font.family"] = "Times New Roman"

#---------------------------------------------------------------------------------

# Functions to read in the raw EDGE galaxy data

# LORES FUNCTION
def ReadGalaxyDataLoRes(galaxy, datatype, downsize='no'):
    
    data = [[] for _ in range(7)]
    # Adding option to read in downsized data if preferred
    if downsize == 'no':
        file_path = f'{galaxy}/{galaxy}_fiducial_{datatype}.data'
    elif downsize == 'yes':
        file_path = f'{galaxy}/{galaxy}_fiducial_downsize_{datatype}.data'
    else:
        raise ValueError("Invalid value for 'downsize'. Choose 'yes' or 'no'.")

    # Assorting the galaxy from the .data files into dictionaries
    with open(file_path, 'r') as file:
        for _ in range(2):
            next(file)
        for line in file:
            values = line.strip().split(',')
            for i in range(7):
                data[i].append(float(values[i]))
                
    print(f"LoRes data for galaxy {galaxy}, {datatype} type, has been read.")
    return data

# HIRES FUNCTION
def ReadGalaxyDataHiRes(galaxy, datatype, downsize='no'):
    
    data = [[] for _ in range(7)]
    # Adding option to read in downsized data if preferred
    if downsize == 'no':
        file_path = f'{galaxy}/{galaxy}_fiducial_hires_{datatype}.data'
    elif downsize == 'yes':
        file_path = f'{galaxy}/{galaxy}_fiducial_hires_downsize_{datatype}.data'
    else:
        raise ValueError("Invalid value for 'downsize'. Choose 'yes' or 'no'.")

    with open(file_path, 'r') as file:
        for _ in range(2):
            next(file)
        for line in file:
            values = line.strip().split(',')
            for i in range(7):
                data[i].append(float(values[i]))
                
        print(f"HiRes data for galaxy {galaxy}, {datatype} type, has been read.")
    return data

#---------------------------------------------------------------------------------
# Functions to read in reduced HiRes data for analysis here

def ReadGobsGbar(galaxies):
    gbar_HiRes = {}
    gobs_HiRes = {}

    for galaxy in galaxies:
        datafile = np.loadtxt(f'{galaxy}/{galaxy}_g_data_HIRES.txt') 

        # Read gbar and gobs data from files
        gbar_data = datafile[:,1]
        gobs_data = datafile[:,2]

        # Store gbar and gobs data in separate dictionaries
        gbar_HiRes[galaxy] = gbar_data
        gobs_HiRes[galaxy] = gobs_data

    return gbar_HiRes, gobs_HiRes

# Defining function to read in galaxy cumulative mass and radius to ensure same profile to LoRes 
# (DM)
def ReadHiResCumMassDM(galaxy):
    data = np.loadtxt(f'{galaxy}/HiRes CumMass Files/{galaxy}_dm_mass.txt')
    radii = data[:,0]
    CumMass = data[:,1]
    
    return CumMass,radii

# (STAR)
def ReadHiResCumMassSTAR(galaxy):
    data = np.loadtxt(f'{galaxy}/{galaxy}_star_mass.txt')
    radii = data[:,0]
    CumMass = data[:,1]
    
    return CumMass,radii

# (GAS)
def ReadHiResCumMassGAS(galaxy):
    data = np.loadtxt(f'{galaxy}/{galaxy}_gas_mass.txt')
    radii = data[:,0]
    CumMass = data[:,1]
    
    return CumMass,radii
#---------------------------------------------------------------------------------

# Function to calculate g baryonic given the galaxy name (e.g. 600) and the dictionary to take data from (LoRes or HiRes)
def CalculateGbaryonic(galaxy, data):
    galaxy_data = data[galaxy]
    stars_data = galaxy_data['star']
    gas_data = galaxy_data['gas']
    
    starsx = np.array(stars_data[0])
    starsy = np.array(stars_data[1])
    starsz = np.array(stars_data[2])
    starmass = np.array(stars_data[6])

    gasx = np.array(gas_data[0])
    gasy = np.array(gas_data[1])
    gasz = np.array(gas_data[2])
    gasmass = np.array(gas_data[6])
    
    G = 6.67e-11
    starsr = np.sqrt(starsx**2 + starsy**2 + starsz**2)
    gasr = np.sqrt(gasx**2 + gasy**2 + gasz**2)
    rSteps = np.logspace(np.log10(100), np.log10(5000), 25)
    CombinedR = np.concatenate((starsr, gasr))
    CombinedMass = np.concatenate((starmass, gasmass))
    
    SortedIndices = np.argsort(CombinedR)
    SortedR = CombinedR[SortedIndices]
    SortedMass = CombinedMass[SortedIndices]
    
    cumulativeMassBaryonic = np.zeros_like(rSteps)
    gbar = np.zeros_like(rSteps)
    for i in range(len(rSteps)):
        particles_within_step = SortedMass[SortedR <= rSteps[i]]
        cumulativeMassBaryonic[i] = np.sum(particles_within_step)
        gbar[i] = (G * cumulativeMassBaryonic[i] * 1.989e30) / (rSteps[i] * 3.086e16)**2

    return gbar

#---------------------------------------------------------------------------------

# Similar function to calculate g observed
def CalculateGobserved(galaxy, data):
    galaxy_data = data[galaxy]
    stars_data = galaxy_data['star']
    gas_data = galaxy_data['gas']
    dm_data = galaxy_data['dm']
    
    starsx = np.array(stars_data[0])
    starsy = np.array(stars_data[1])
    starsz = np.array(stars_data[2])
    starmass = np.array(stars_data[6])

    gasx = np.array(gas_data[0])
    gasy = np.array(gas_data[1])
    gasz = np.array(gas_data[2])
    gasmass = np.array(gas_data[6])
    
    DMx = np.array(dm_data[0])
    DMy = np.array(dm_data[1])
    DMz = np.array(dm_data[2])
    DMmass = np.array(dm_data[6])
    
    G = 6.67e-11
    starsr = np.sqrt(starsx**2 + starsy**2 + starsz**2)
    gasr = np.sqrt(gasx**2 + gasy**2 + gasz**2)
    DMr = np.sqrt(DMx**2 + DMy**2 + DMz**2)
    rSteps = np.logspace(np.log10(100), np.log10(5000), 25)

    CombinedR = np.concatenate((starsr, gasr, DMr))
    CombinedMass = np.concatenate((starmass, gasmass, DMmass))
    
    SortedIndices = np.argsort(CombinedR)
    SortedR = CombinedR[SortedIndices]
    SortedMass = CombinedMass[SortedIndices]
    
    cumulativeMass = np.zeros_like(rSteps)
    gobs = np.zeros_like(rSteps)
    for i in range(len(rSteps)):
        particles_within_step = SortedMass[SortedR <= rSteps[i]]
        cumulativeMass[i] = np.sum(particles_within_step) #Use indexing instead
        gobs[i] = (G * cumulativeMass[i] * 1.989e30) / (rSteps[i] * 3.086e16)**2

    return gobs

#---------------------------------------------------------------------------------

# Function that calculates mass profile within galaxy
def CalculateMassProfile(galaxy, data, particletype):
    galaxy_data = data[galaxy]

    particle_data = galaxy_data[particletype]
    x = np.array(particle_data[0])
    y = np.array(particle_data[1])
    z = np.array(particle_data[2])
    mass = np.array(particle_data[6])

    r = np.sqrt(x**2 + y**2 + z**2)
    sorted_indices = np.argsort(r)
    sorted_r = r[sorted_indices]
    sorted_mass = mass[sorted_indices]

    cumulative_mass = np.cumsum(sorted_mass)

    return cumulative_mass, sorted_r

#---------------------------------------------------------------------------------

# Latest dwarf galaxy experimental data (see Justin Read et al, MNRAS, 484 (2019))
# List of galaxies
Exp_galaxies = ['AntliaB','Carina','Draco','Eri2','Fornax','Grus1','HydraII','LeoT','LI','LII','Sextans','UMi']

def read_experimental_data(galaxies):
    ExperimentalData = {}

    for galaxy in Exp_galaxies:
        file_path = f'Experimental Data/GbarGobsExpData/{galaxy}_gbar_gobs_data.txt'
        with open(file_path, 'r') as file:
            data = []
            for line in file:
                radius, gbar, gobs, gobs_low, gobs_high = map(float, line.split())
                data.append((radius, gbar, gobs, gobs_low, gobs_high))
            ExperimentalData[galaxy] = data

    return ExperimentalData

ExperimentalData = read_experimental_data(Exp_galaxies)

#---------------------------------------------------------------------------------

# Reading in data for both LoRes and HiRes simulations
galaxies = [600, 605, 624, 1445, 1459]
LoResData = {}
HiResData = {}

# Read data for all galaxies
for galaxy in galaxies:
    LoResData[galaxy] = {}
    HiResData[galaxy] = {}
    for datatype in ['star', 'gas', 'dm']:
        LoResData[galaxy][datatype] = ReadGalaxyDataLoRes(galaxy, datatype, 'no')
        HiResData[galaxy][datatype] = ReadGalaxyDataHiRes(galaxy, datatype, 'yes')

# #---------------------------------------------------------------------------------

# Function to find total mass of substance of galaxy
def calculate_total_mass(particle_data):
    total_mass = sum(particle_data)
    return total_mass
 

# Full-Size HiRes masses supplied by Justin Read originally used to correct the total masses of the downsized files
HiResFullSizeMass = {
    624: {'dm': 10**9.414926244657014, 'star': 10**6.034727296506947, 'gas': 10**7.0624389720272305},
    600: {'dm': 10**9.506971788041424, 'star': 10**5.993043811198278, 'gas': 10**6.926655490229328},
    605: {'dm': 10**9.502840380683011, 'star': 10**6.285702674992516, 'gas': 10**6.922411072805887},
    1445: {'dm': 10**9.117011822623084, 'star': 10**5.129770291886941, 'gas': 10**6.373640127139425},
    1459: {'dm': 10**9.152991186834008, 'star': 10**5.576396672565849, 'gas': 10**6.465697114988802}  # No gas data for galaxy 1459 in downsized files
}


# Loop through all galaxies to correct for total mass in downsized files
for galaxy in galaxies:
    # Loop through particle types (dm and gas)
    for particle_type in ['dm', 'gas']:
        # Calculate total mass of the downsized HiRes data
        total_mass_downsized_hires = calculate_total_mass(HiResData[galaxy][particle_type][6])

        # Determine the appropriate scaling factor based on particle type and galaxy
        if galaxy != 1459:
            if particle_type == 'dm':
                scaling_factor = HiResFullSizeMass[galaxy]['dm'] / total_mass_downsized_hires
            elif particle_type == 'gas':
                scaling_factor = HiResFullSizeMass[galaxy]['gas'] / total_mass_downsized_hires
        else:
            # Skip adjustment for gas content in galaxy 1459
            if particle_type == 'dm':
                scaling_factor = HiResFullSizeMass[galaxy]['dm'] / total_mass_downsized_hires
            else:
                scaling_factor = 1.0

        # Adjust individual particle masses in the downsized HiRes data using the scaling factor
        particle_masses_downsized_hires_adjusted = [mass * scaling_factor for mass in HiResData[galaxy][particle_type][6]]

        # Update HiResData with adjusted particle masses
        HiResData[galaxy][particle_type][6] = particle_masses_downsized_hires_adjusted
        

# These calculations are no longer necessary and are bypassed due to the calculation of the full-size HiRes files seperately (see GITHUB page)
#---------------------------------------------------------------------------------

# Plotting X vs. Y for stars data in galaxies
fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex='col', sharey='row')

for i, galaxy in enumerate(galaxies):
    # Plot LoResData
    starsx_lores = LoResData[galaxy]['star'][0]
    starsy_lores = LoResData[galaxy]['star'][1]
    H_lores, xedges_lores, yedges_lores = np.histogram2d(starsx_lores, starsy_lores, bins=50)
    axes[0, i].tick_params(axis='both', which='major', labelsize=20)
    axes[0, i].set_title(f'{galaxy}', fontsize=20)
    axes[0, i].contourf(xedges_lores[:-1], yedges_lores[:-1], np.log10(H_lores.T + 1), cmap='Blues', alpha=0.7, levels=np.linspace(0, np.log10(H_lores.max() + 1), 20))
    
    # Plot HiResData
    starsx_hires = HiResData[galaxy]['star'][0]
    starsy_hires = HiResData[galaxy]['star'][1]
    H_hires, xedges_hires, yedges_hires = np.histogram2d(starsx_hires, starsy_hires, bins=50)
    axes[1, i].tick_params(axis='both', which='major', labelsize=20)
    axes[1, i].contourf(xedges_hires[:-1], yedges_hires[:-1], np.log10(H_hires.T + 1), cmap='Reds', alpha=0.7, levels=np.linspace(0, np.log10(H_hires.max() + 1), 20))

    # Set x-axis ticks for the bottom row
    if i == 4:
        axes[1, i].set_xticks([-2000, -1000, 0, 1000, 2000])
        axes[1, i].set_xticklabels(['-2', '-1', '0', '1', '2'], fontsize=18)
    else:
        axes[1, i].set_xticks([-2000, -1000, 0, 1000])
        axes[1, i].set_xticklabels(['-2', '-1', '0', '1'], fontsize=18)
    
    # Set y-axis ticks and labels for the bottom row
    if i == 0:
        axes[1, i].set_yticks([-2000, -1000, 0, 1000, 2000])
        axes[1, i].set_yticklabels(['-2', '-1', '0', '1', '2'], fontsize=18)
        axes[1, i].set_ylabel('Y (kpc)', fontsize=18)
    else:
        axes[1, i].set_yticks([-2000, -1000, 0, 1000])
        axes[1, i].set_yticklabels(['-2', '-1', '0', '1'], fontsize=18)
    
    # Set y-axis ticks and labels for the top row
    if i == 0:
        axes[0, i].set_yticks([-2000, -1000, 0, 1000, 2000])
        axes[0, i].set_yticklabels(['-2', '-1', '0', '1', '2'], fontsize=18)
        axes[0, i].set_ylabel('Y (kpc)', fontsize=18)
    else:
        axes[0, i].set_yticks([-2000, -1000, 0, 1000,2000])
        axes[0, i].set_yticklabels(['-2', '-1', '0', '1','2'], fontsize=18)
    
# Share x-axis labels
for ax in axes[1, :]:
    ax.set_xlabel('X (kpc)', fontsize=18)

# Set x-axis limits for bottom row
for ax in axes[1, :]:
    ax.set_xlim(-2000, 2000)

# Set y-axis limits for left column
for ax in axes[:, 0]:
    ax.set_ylim(-2000, 2000)

plt.subplots_adjust(hspace=0, wspace=0)  # Remove horizontal and vertical space between subplots
plt.savefig('EDGE_XY_Stars.png', dpi=610)
plt.show()

#---------------------------------------------------------------------------------

fig, axs = plt.subplots(len(galaxies), 1, figsize=(8, 12), sharex=True, sharey=True)

# Loop over all 5 galaxies
for i, galaxy in enumerate(galaxies):
    # Calculate mass profile
    cumulative_mass, sorted_r = CalculateMassProfile(galaxy, LoResData,'gas')
    HiRes_cum_mass, HiRes_sorted_r = ReadHiResCumMassSTAR(galaxy)
    
    # Plot
    axs[i].plot(sorted_r, cumulative_mass, color='slategrey', linestyle='-', label='LoRes Data', alpha=0.8)
    axs[i].plot(HiRes_sorted_r, HiRes_cum_mass, color='mediumturquoise', linestyle='--', label='HiRes Data', alpha=0.8)
    axs[i].set_title(f'Galaxy {galaxy}', fontsize=12)
    axs[i].set_xlabel('Radius (pc)', fontsize=10)
    axs[i].set_ylabel('Cumulative Mass ($M_{\odot}$)', fontsize=10)
    axs[i].set_xlim(0,30000)
    axs[i].set_xscale('log')
    axs[i].tick_params(axis='both', which='major', labelsize=8)
    axs[i].legend(fontsize=8)

plt.tight_layout()

plt.subplots_adjust(hspace=0.3)
plt.savefig('GASmassprofilesTROUBLESHOOT.png', dpi=610) 
plt.show()

#---------------------------------------------------------------------------------

# Calculating and plotting gbar and gobs for all galaxies
gbar_LoRes = {}
gobs_LoRes = {}

for galaxy in galaxies:
    gbar_LoRes[galaxy] = CalculateGbaryonic(galaxy, LoResData)
    gobs_LoRes[galaxy] = CalculateGobserved(galaxy, LoResData)

# READING IN NEW DATA FROM EXTERNALLY CALCULATED GBAR AND GOBS
gbar_HiRes, gobs_HiRes = ReadGobsGbar(galaxies)

# #---------------------------------------------------------------------------------

# USE THIS IF USING DOWNSIZED DATA!!!
# # Calculate Gbaryonic and Gobserved for HiRes data (OLD VERSION)
# gbar_HiRes = {}
# gobs_HiRes = {}
# for galaxy in galaxies:
#     gbar_HiRes[galaxy] = CalculateGbaryonic(galaxy, HiResData)
#     gobs_HiRes[galaxy] = CalculateGobserved(galaxy, HiResData)

#-----------------------------------------------------------------------------------

# Calculations below show the RAR from the Lelli + McGaugh paper about RAR (2017), used to compare with our simulation RAR data

# Defining functions for RAR Calculations for plotting analytical relations later
def LTGRARCalc(x,gcross):
    y = (x / (1 - np.exp(-np.sqrt(x / (gcross)))))
    return y

# Completing all calculations for RAR analytical curves and their associated errors
# Defining parameters for RAR curves (Lelli 2017)
RARxData = np.logspace(-16, -10, 100)
gcrossLTG = 1.2e-10
gcrossLTGerr = 0.02e-10

ghat = 9.2e-12
ghaterr = 0.2e-12 

# Generate samples using gaussian, to calculate average uncertainty on RAR
gcrossLTGSamples = np.random.normal(gcrossLTG, gcrossLTGerr, 100)
ghatSamples = np.random.normal(ghat, ghaterr, 100)

LTGRARSample = np.array([LTGRARCalc(RARxData, gcrossLTG_sample) for gcrossLTG_sample in gcrossLTGSamples])

# Calculate mean and standard deviation of RAR random samples
LTGRARmean = np.mean(LTGRARSample, axis=0)
LTGRARstd = np.std(LTGRARSample, axis=0)
LTGRARUpperBound = LTGRARmean + 1*LTGRARstd #Change multiplier of std for 68% criterion (1*), 99.7% criterion (3*)
LTGRARLowerBound = LTGRARmean - 1*LTGRARstd #Change multiplier of std for 68% criterion (1*), 99.7% criterion (3*)

# Calculate true values
LTGRARTrue = LTGRARCalc(RARxData, gcrossLTG)

# # Plot RAR if needed:
# fig, axs = plt.subplots(1, 2, figsize=(15, 5))

# # Plot individual samples
# axs[0].plot(RARxData, LTGRARSample.T, color='red', alpha=0.1, label='LTGRAR Sample')
# axs[0].set_xscale('log')
# axs[0].set_yscale('log')
# axs[0].set_xlabel('gbar (m/s^2)')
# axs[0].set_ylabel('gobs (m/s^2)')
# axs[0].set_title('Individual Samples')

# # Create legend handles for only the desired lines
# blue_line = plt.Line2D([], [], color='red', label='LTGRAR Sample')
# # Display legend with only the specified lines
# axs[0].legend(handles=[blue_line, green_line])

# # Plot histograms
# axs[1].hist(LTGRARSample.flatten(), bins=20, color='red', alpha=0.7, label='LTGRAR Sample')
# axs[1].set_xlabel('Value')
# axs[1].set_ylabel('Frequency')
# axs[1].set_title('Histograms')
# axs[1].legend()

# plt.tight_layout()
# plt.savefig('ErroronRAR.png', dpi=610)
# plt.show()

#---------------------------------------------------------------------------------

# Plotting the real dwarf galaxy data from Read et al (2019):
custom_palette = ['mediumseagreen', 'mediumturquoise', 'cornflowerblue', 'palegreen', 'slategrey', 'lightcoral', 'goldenrod', 'orchid', 'deepskyblue', 'olivedrab']

# # Create a new figure
plt.figure(figsize=(8, 8))  # Changed figsize to make square plot

# Plot experimental data for each galaxy using custom colors
for i, galaxy in enumerate(Exp_galaxies):
    data = ExperimentalData[galaxy]
    radius, gbar, gobs, gobs_low, gobs_high = zip(*data)
    color = custom_palette[i % len(custom_palette)]  # Ensure cycling through colors if more galaxies than colors
    plt.scatter(gbar, gobs, color=color, s=4)
    plt.fill_between(gbar, gobs_low, gobs_high, alpha=0.4, color=color, label=f'{galaxy}')

# Plot 1:1 line
x = np.logspace(-16, -10, 100)
plt.plot(x, x, 'k--', label='1:1 line')

# Plot RAR analytical relationship if wanted
# plt.plot(RARxData, LTGRARTrue, label='MOND (McGaugh 2016)', color='red', linestyle=':')
# plt.fill_between(RARxData, LTGRARLowerBound, LTGRARUpperBound, color='pink', alpha=0.4)

# Set labels and scale
plt.xlabel('$g_{\mathrm{bar}}$ (m/s$^2$)', fontsize=16)
plt.xscale('log')
plt.ylabel('$g_{\mathrm{obs}}$ (m/s$^2$)', fontsize=16)
plt.yscale('log')
plt.xlim(1e-15, 1e-11)
plt.ylim(1e-14, 1e-10)  # Changed to same range as x-axis
plt.tick_params(axis='both', which='major', labelsize=16)

# Add legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=14)

# Set aspect ratio to be equal
plt.gca().set_aspect('equal')

# Title and tight layout
# plt.title('Experimental Data: gbar vs gobs with 1$\sigma$ uncertainties')
plt.tight_layout()

# Save or show the plot
plt.savefig('TESTexperimental_gbar_vs_gobs.png', dpi=610)
plt.show()

#---------------------------------------------------------------------------------

# Plot individual plots for the 5 galaxies. Each plot shows HiRes and LoRes EDGE data

markers = ['o', 's', '^', 'D', 'P']
colors = ['mediumseagreen', 'mediumturquoise', 'cornflowerblue', 'palegreen', 'slategrey']
colors2 = ['springgreen','paleturquoise','royalblue','forestgreen','gainsboro']

# Plot LoRes and HiRes data for each galaxy
for i, galaxy in enumerate(galaxies):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Plot LoRes data
    ax.scatter(gbar_LoRes[galaxy], gobs_LoRes[galaxy], label=f'LoRes ($M_{{\mathrm{{DM}}}}$ = 1112 $M_{{\odot}}$)', marker='o', color=colors[i], edgecolor='black')
    
    # Plot HiRes data
    ax.scatter(gbar_HiRes[galaxy], gobs_HiRes[galaxy], label=f'HiRes ($M_{{\mathrm{{DM}}}}$ = 117 $M_{{\odot}}$)', marker='s', color=colors2[i], edgecolor='black')
    
    # Plot 1:1 line
    x = np.logspace(-16, -10, 100)
    y = np.logspace(-16, -10, 100)
    ax.plot(x, y, 'k--', label='1:1 line')
    
    # Plot RAR analytical relationship
    ax.plot(RARxData, LTGRARTrue, label='MOND (McGaugh 2016)', color='red', linestyle=':')
    ax.fill_between(RARxData, LTGRARLowerBound, LTGRARUpperBound, color='pink', alpha=0.4)
    
    # Set labels and scales
    # ax.set_title(f'Galaxy {galaxy}')
    ax.set_xlabel('$g_{\mathrm{bar}}$ (m/s$^2$)',fontsize=17)
    ax.set_xscale('log')
    ax.set_ylabel('$g_{\mathrm{obs}}$ (m/s$^2$)',fontsize=17)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    ax.text(0.05, 0.95, '100 pc < r < 5000 pc', transform=ax.transAxes, fontsize=17, verticalalignment='top')
    ax.legend(fontsize=14)
        
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(f'Galaxy_{galaxy}_gbargobs.png', dpi=610)
    plt.show()

#---------------------------------------------------------------------------------

# Plotting all LoRes and HiRes EDGE data on respective plots to gauge overall distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# # Define markers and colors for galaxies
markers = ['o', 's', '^', 'D', 'P']
colors = ['mediumseagreen', 'mediumturquoise', 'cornflowerblue', 'palegreen', 'slategrey']
colors2 = ['green','paleturquoise','blue','limegreen','gainsboro']

# Plot LoRes data
axes[0].set_title('LoRes Data')
for i, galaxy in enumerate(galaxies):
    axes[0].scatter(gbar_LoRes[galaxy], gobs_LoRes[galaxy], label=f'Galaxy {galaxy}', marker=markers[i], color=colors[i], edgecolor='black')
    # axes[0].plot(gbar_LoRes[galaxy], gobs_LoRes[galaxy], linewidth=2, color=colors[i])

# Plot 1:1 line for LoRes
x = np.logspace(-16, -10, 100)
y = np.logspace(-16, -10, 100)
axes[0].plot(x, y, 'k--', label='1:1 line')

# Plot RAR analytical relationship for LoRes
axes[0].plot(RARxData, LTGRARTrue, label='MOND (McGaugh 2016)', color='red', linestyle=':')
axes[0].fill_between(RARxData, LTGRARLowerBound, LTGRARUpperBound, color='pink', alpha=0.4)

# Plotting experimental data on plot
for galaxy in Exp_galaxies:
    data = ExperimentalData[galaxy]
    radius, gbar, gobs, gobs_low, gobs_high = zip(*data)
    axes[0].fill_between(gbar, gobs_low, gobs_high, color='gray', alpha=0.4)
axes[0].scatter([], [], color='gray', label='Experimental Data')

axes[0].set_xlabel('$g_{\mathrm{bar}}$ (m/s$^2$)',fontsize=15)
axes[0].set_xscale('log')
axes[0].set_ylabel('$g_{\mathrm{obs}}$ (m/s$^2$)',fontsize=15)
axes[0].set_yscale('log')
axes[0].tick_params(axis='both', which='major', labelsize=15)

axes[0].text(0.05, 0.95, '100 pc < r < 5000 pc', transform=axes[0].transAxes, fontsize=14, verticalalignment='top')
axes[0].text(0.05, 0.90, '$M_{\mathrm{DM}}$ = 1112 $M_{\odot}$', transform=axes[0].transAxes, fontsize=14, verticalalignment='top')
axes[0].legend(fontsize=13)

# Plot HiRes data
axes[1].set_title('HiRes Data')
for i, galaxy in enumerate(galaxies):
    axes[1].scatter(gbar_HiRes[galaxy], gobs_HiRes[galaxy], label=f'Galaxy {galaxy}', marker=markers[i], color=colors[i], edgecolor='black')
    # axes[1].plot(gbar_HiRes[galaxy], gobs_HiRes[galaxy], linewidth=2, color=colors[i])

# Plot 1:1 line for HiRes
axes[1].plot(x, y, 'k--', label='1:1 line')

# Plot RAR analytical relationship for HiRes
axes[1].plot(RARxData, LTGRARTrue, label='MOND (McGaugh 2016)', color='red', linestyle=':')
axes[1].fill_between(RARxData, LTGRARLowerBound, LTGRARUpperBound, color='pink', alpha=0.4)

for galaxy in Exp_galaxies:
    data = ExperimentalData[galaxy]
    radius, gbar, gobs, gobs_low, gobs_high = zip(*data)
    axes[1].fill_between(gbar, gobs_low, gobs_high, color='gray', alpha=0.4)
axes[1].scatter([], [], color='gray', label='Experimental Data')

axes[1].set_xlabel('$g_{\mathrm{bar}}$ (m/s$^2$)',fontsize=15)
axes[1].set_xscale('log')
axes[1].set_ylabel('$g_{\mathrm{obs}}$ (m/s$^2$)',fontsize=15)
axes[1].set_yscale('log')
axes[1].tick_params(axis='both', which='major', labelsize=15)
axes[1].text(0.05, 0.95, '100 pc < r < 5000 pc', transform=axes[1].transAxes, fontsize=14, verticalalignment='top')
axes[1].text(0.05, 0.90, '$M_{\mathrm{DM}}$ = 117 $M_{\odot}$', transform=axes[1].transAxes, fontsize=14, verticalalignment='top')
axes[1].legend(fontsize=13)

axes[0].set_aspect('equal', adjustable='box')
axes[1].set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig('gbargobsWITHexperimental.png', dpi=610)
plt.show()

#---------------------------------------------------------------------------------

# Plotting Residuals of all EDGE data and Experimental data, from the MONDian RAR (Lelli + McGaugh)
plt.figure(figsize=(12, 6))

# Plot HiRes data with residuals
for galaxy in galaxies:
    RAR_residuals = np.array(gobs_HiRes[galaxy]) - np.interp(gbar_HiRes[galaxy], RARxData, LTGRARTrue)
    plt.scatter(gbar_HiRes[galaxy], RAR_residuals, marker='^', color='coral', label='_nolegend_',alpha=0.5,edgecolors='black',s=50)

# Plot LoRes data with residuals
for galaxy in galaxies:
    RAR_residuals = np.array(gobs_LoRes[galaxy]) - np.interp(gbar_LoRes[galaxy], RARxData, LTGRARTrue)
    plt.scatter(gbar_LoRes[galaxy], RAR_residuals, marker='^', color='mediumturquoise', label='_nolegend_',alpha=0.5,edgecolors='black',s=50)

# Plot experimental data points with residuals
for galaxy in Exp_galaxies:
    data = ExperimentalData[galaxy]
    radius, gbar, gobs, _, _ = zip(*data)
    RAR_residuals = np.array(gobs) - np.interp(gbar, RARxData, LTGRARTrue)
    plt.scatter(gbar, RAR_residuals, color='gray', alpha=0.3,edgecolors='black',marker='o')

plt.scatter([], [], color='mediumturquoise', marker='^', label='EDGE LoRes Data',edgecolors='black')
plt.scatter([], [], color='coral', marker='^', label='EDGE HiRes Data',edgecolors='black')

plt.scatter([], [], color='gray', label='Experimental Data (Read et al)')
plt.axhline(y=0, color='red', linestyle='--', linewidth=2, label='MOND (McGaugh 2016)')
plt.axhline(y=0.5e-11, color='black', linestyle=':', linewidth=2)
plt.axhline(y=-0.5e-11, color='black', linestyle=':', linewidth=2)

plt.xlabel('$g_{\mathrm{bar}}$ (m/s$^2$)', fontsize=16)
plt.ylabel('Residuals', fontsize=16)
plt.xscale('log')
plt.ylim(-2.5e-11, 2.5e-11)
plt.xlim(1e-15, 1e-11)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=16)

plt.tight_layout()
plt.savefig('Residuals_All_Data_Plot_TEST.png', dpi=610)
plt.show()

#---------------------------------------------------------------------------------

# Saving Gbar and Gobs data to textfiles for further analysis. To be published in M. Julio et al. (2024)
# Loop through all galaxies
for galaxy in galaxies:
    # Saving gbar and gobs data for LoRes
    with open(f'EDGEGalaxy{galaxy}LoResGbarGobs.txt', 'w') as file:
        file.write("Galaxy | Radius (pc) | gobs (ms^-2) | gbar (ms^-2)\n")
        gbar_data = gbar_LoRes[galaxy]
        gobs_data = gobs_LoRes[galaxy]
        radius_data = np.logspace(np.log10(100), np.log10(5000), 25)  # Assuming radius data is the same for all galaxies
        for radius, gobs_val, gbar_val in zip(radius_data, gobs_data, gbar_data):
            file.write(f"{galaxy}  {radius}  {gobs_val}  {gbar_val}\n")
    
    # Saving gbar and gobs data for HiRes
    with open(f'EDGEGalaxy{galaxy}HiResGbarGobs.txt', 'w') as file:
        file.write("Galaxy | Radius (pc) | gobs (ms^-2) | gbar (ms^-2)\n")
        gbar_data = gbar_HiRes[galaxy]
        gobs_data = gobs_HiRes[galaxy]
        radius_data = np.logspace(np.log10(100), np.log10(5000), 25)  # Assuming radius data is the same for all galaxies
        for radius, gobs_val, gbar_val in zip(radius_data, gobs_data, gbar_data):
            file.write(f"{galaxy}  {radius}  {gobs_val}  {gbar_val}\n")

#---------------------------------------------------------------------------------
