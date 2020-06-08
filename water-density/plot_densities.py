
import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec


# Set font sizes for the graph
small_font = 10
medium_font = 12
large_font = 14

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=small_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=small_font)


def read_csv(filename):
    """
    Read in density data from a CSV file

    Parameters
    ----------
    filename : str
        Name of CSV file

    Returns
    -------
    times : list
        Time values for each density
    means : numpy.ndarray
        Mean density at each point in time
    errors : numpy.ndarray
        Std. error in mean density at each point in time
    densities : list
        All density values observed
    """
    # Read lines from file
    with open(filename, 'r') as f:
        lines = f.readlines()

    times = []
    means = []
    errors = []
    densities = []

    for line in lines[1:]:
        # Read time
        t = float(line.split(',')[0])
        times.append(t)
        # Read the densities
        line_densities = [float(x) for x in line.split(',')[1:]]
        means.append(np.mean(line_densities))
        errors.append(np.std(line_densities) / np.sqrt(len(line_densities)))
        # Store each individual density
        for d in line_densities:
            densities.append(d)

    # Convert means and errors to numpy arrays
    means = np.array(means)
    errors = np.array(errors)

    return times, means, errors, densities


parser = argparse.ArgumentParser()
parser.add_argument('--npt', default=None, help='CSV containing densities from NPT simulation')
parser.add_argument('--uvt', default=None, help='CSV containing densities from GCMC/MD simulation')
parser.add_argument('--xlim', default=None, type=float, help='Max value for the x axis')
parser.add_argument('--ylims', default=None, type=float, nargs='+', help='Min and max values for the y axis')
parser.add_argument('-o', '--output', default='fig.pdf', help='Name of the output file')
args = parser.parse_args()

# Read NPT densities
if args.npt is not None:
    # Read file
    npt_times, npt_means, npt_errors, npt_densities = read_csv(args.npt)
    # Print out overall mean and std. deviation
    npt_mean = np.mean(npt_densities)
    npt_std_dev = np.std(npt_densities)
    print('NPT density = {:.6f} g/mL (SD = {:.6f} g/mL)'.format(np.round(npt_mean, 6), np.round(npt_std_dev, 6)))

# Read in GCMC/MD densities
if args.uvt is not None:
    # Read file
    uvt_times, uvt_means, uvt_errors, uvt_densities = read_csv(args.uvt)
    # Print out overall mean and std. deviation
    uvt_mean = np.mean(uvt_densities)
    uvt_std_dev = np.std(uvt_densities)
    print('uVT density = {:.6f} g/mL (SD = {:.6f} g/mL)'.format(np.round(uvt_mean, 6), np.round(uvt_std_dev, 6)))

# Configure graph
fig = plt.figure(figsize=(14,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Plot NPT data
if args.npt is not None:
    # Plot density vs time
    ax1.plot(npt_times, npt_means, color='blue', alpha=1.0, label=r'NPT')
    # Shade the area within 1 standard error of the mean
    ax1.fill_between(npt_times, npt_means-npt_errors, npt_means+npt_errors, lw=0, color='blue', alpha=0.3)

    # Plot density histogram
    ax2.hist(npt_densities, bins='auto', color='blue', density=True, orientation='horizontal', histtype='stepfilled', alpha=0.3, lw=0)
    ax2.hist(npt_densities, bins='auto', color='blue', density=True, orientation='horizontal', histtype='step', alpha=1.0, lw=1.5)

# Plot uVT data
if args.uvt is not None:
    # Plot density vs time
    ax1.plot(uvt_times, uvt_means, color='red', alpha=1.0, label=r'$\mu$VT')
    # Shade the area within 1 standard error of the mean
    ax1.fill_between(uvt_times, uvt_means-uvt_errors, uvt_means+uvt_errors, lw=0, color='red', alpha=0.3)
    
    # Plot density histogram
    ax2.hist(uvt_densities, bins='auto', color='red', density=True, orientation='horizontal', histtype='stepfilled', alpha=0.3, lw=0)
    ax2.hist(uvt_densities, bins='auto', color='red', density=True, orientation='horizontal', histtype='step', alpha=1.0, lw=1.5)

# Set upper x-limit
if args.xlim is not None:
    ax1.set_xlim(left=0, right=args.xlim)
else:
    ax1.set_xlim(left=0)

# Set y-limits
if args.ylims is not None:
    ax1.set_ylim(args.ylims[0], args.ylims[1])
    ax2.set_ylim(args.ylims[0], args.ylims[1])

ax2.set_xticks([])
ax2.set_yticks([])
ax1.set_xlabel('Time / ns')
ax1.set_ylabel(r'Density / g mL$^{-1}$')
ax1.legend()
ax1.axis('on')
ax2.axis('off')
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.savefig(args.output)
plt.show()

