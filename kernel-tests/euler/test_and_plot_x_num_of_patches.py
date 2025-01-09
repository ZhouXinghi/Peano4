import os
import subprocess
import re
import matplotlib.pyplot as plt
import json

def run_test(ps_value, d_value):
    # Clean previous builds
    subprocess.run("make clean && make distclean", shell=True, check=True)
    
    # Run the Python script with the specified -ps and -d values
    # This generates the C++ files required for benchmarking
    command = f"python kernel-benchmarks-fv-rusanov.py -d {d_value} -ps {ps_value} -gpu --iterations 10 --ite_per_transfer 1"
    subprocess.run(command, shell=True, check=True)

    # Compile the generated C++ files using the cmake.sh script
    subprocess.run("sh cmake.sh", shell=True, check=True)

    # Run the compiled executable and capture the output for further analysis
    result = subprocess.run("build/kernelTestMain", shell=True, check=True, capture_output=True, text=True)
    return result.stdout

def parse_results(output):
    # Use regex to extract timing information from the output of the executable
    pattern = r"========Patches = (\d+)\n.*?GPUBaseline: (\d+\.\d+)s.*?GPUPacked: (\d+\.\d+)s.*?"
    matches = re.findall(pattern, output, re.DOTALL)
    
    # Store the extracted data in a list of dictionaries for easy access
    results = []
    for match in matches:
        patches, gpu_baseline, gpu_packed = match
        results.append({
            'patches': int(patches),
            'gpu_baseline': float(gpu_baseline),
            'gpu_packed': float(gpu_packed),
        })
    return results

def plot_results(results, ps_values, d_value, output_folder):
    # Plot the results for each patch size value
    # output_folder = "results"
    os.makedirs(output_folder, exist_ok=True)
    
    for ps_value, result_set in zip(ps_values, results):
        # Extract data for plotting
        patches = [r['patches'] for r in result_set]
        # cpu_times = [r['cpu'] for r in result_set]
        gpu_baseline_times = [r['gpu_baseline'] for r in result_set]
        gpu_packed_times = [r['gpu_packed'] for r in result_set]

        # Create a plot for the current patch size
        plt.figure(figsize=(10, 6))
        # plt.plot(patches, cpu_times, label='CPU', marker='o')
        plt.plot(patches, gpu_baseline_times, label='GPU Baseline', marker='o')
        plt.plot(patches, gpu_packed_times, label='GPU Packed', marker='o')

        # Add labels, title, legend, and grid to the plot
        plt.xlabel('Number of Patches')
        plt.ylabel('Time (s)')
        plt.title(f'Performance Comparison for Patch Size {ps_value}, Dimension {d_value}')
        plt.legend()
        plt.grid(True)
        # Save the plot as a PNG file in the output folder
        plt.savefig(os.path.join(output_folder, f'results_ps_{ps_value}_d_{d_value}.png'))
        # Do not display the plot


def main():
    d_values = [2]
    results_folder = "results-thesis"
    os.makedirs(results_folder, exist_ok=True)

    for d_value in d_values:
        if d_value == 2:
            ps_values = [4, 8, 12, 16, 20, 24, 28, 32, 36]
        else:
            ps_values = []

        results_file = os.path.join(results_folder, f"test_results_d_{d_value}.json")

        # Check if the results file already exists
        if os.path.exists(results_file):
            # Load results from the existing file
            with open(results_file, 'r') as file:
                all_results = json.load(file)
            print(f"Loaded existing results from file for dimension {d_value}.")
        else:
            # Run tests for each patch size and collect the results
            all_results = []
            for ps in ps_values:
                print(f"Running test for patch size {ps}, dimension {d_value}...")
                output = run_test(ps, d_value)
                results = parse_results(output)
                all_results.append(results)

            # Save the results to a file in the results folder
            with open(results_file, 'w') as file:
                json.dump(all_results, file)
            print(f"Test results saved to file for dimension {d_value}.")

        # Plot the results for all patch sizes for the current dimension value
        plot_results(all_results, ps_values, d_value, results_folder)

if __name__ == "__main__":
    main()
