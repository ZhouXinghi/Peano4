import os
import subprocess
import re
import matplotlib.pyplot as plt
import json

def run_test(ps_value, d_value, ite_per_transfer):
    # Clean previous builds
    subprocess.run("make clean && make distclean", shell=True, check=False)
    # Run the Python script with the specified -ps and -d values, and ite_per_transfer
    command = f"python kernel-benchmarks-fv-rusanov.py -d {d_value} -ps {ps_value} -p 2048 -gpu --iterations 100 --ite_per_transfer {ite_per_transfer}"
    subprocess.run(command, shell=True, check=True)

    # Compile the generated C++ files using the cmake.sh script
    subprocess.run("sh cmake.sh", shell=True, check=True)

    # Run the compiled executable and capture the output for further analysis
    result = subprocess.run("build/kernelTestMain", shell=True, check=True, capture_output=True, text=True)
    return result.stdout

def parse_results(output):
    # Use regex to extract timing information from the output of the executable
    # pattern = r"========Patches = (\d+)\n.*?CPU: (\d+\.\d+)s.*?timeStepWithRusanovPatchwiseHeapStateless: (\d+\.\d+)s.*?PatchWiseGPUPacked: (\d+\.\d+)s"
    pattern = r"========Patches = (\d+)\n.*?GPUBaseline: (\d+\.\d+)s.*?GPUPacked: (\d+\.\d+)s"
    matches = re.findall(pattern, output, re.DOTALL)
    
    # Store the extracted data in a list of dictionaries for easy access
    results = []
    for match in matches:
        patches, gpu_baseline, gpu_packed = match
        results.append({
            'patches': int(patches),
            # 'cpu': float(cpu),
            'gpu_baseline': float(gpu_baseline),
            'gpu_packed': float(gpu_packed)
        })
    return results

def plot_results(all_results, ps_values, ite_per_transfer_values, d_value, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    
    for idx, ite_per_transfer in enumerate(ite_per_transfer_values):
        # Create lists to hold time data for each patch size at the current ite_per_transfer
        cpu_times = []
        gpu_baseline_times = []
        gpu_packed_times = []

        # Gather time data for each patch size at the current ite_per_transfer
        for ps_results in all_results:
            cpu_times.append(ps_results[idx].get('cpu', 0))
            gpu_baseline_times.append(ps_results[idx].get('gpu_baseline', 0))
            gpu_packed_times.append(ps_results[idx].get('gpu_packed', 0))

        # Create a plot for the current ite_per_transfer value
        plt.figure(figsize=(12, 7))
        # plt.plot(ps_values, cpu_times, label='CPU', marker='o')
        plt.plot(ps_values, gpu_baseline_times, label='GPU Baseline', marker='o')
        plt.plot(ps_values, gpu_packed_times, label='GPU Packed', marker='o')

        # Add labels, title, legend, and grid to the plot
        plt.xlabel('Patch Size')
        plt.ylabel('Time (s)')
        plt.title(f'Performance Comparison for ite_per_transfer {ite_per_transfer}, Dimension {d_value}')
        plt.legend()
        plt.grid(True)
        
        # Save the plot as a PNG file in the output folder
        plt.savefig(os.path.join(output_folder, f'results_d_{d_value}_ite_{ite_per_transfer}.png'))
        plt.close()

def main():
    d_values = [2]
    results_folder = "results-hpc-ext-thesis"
    os.makedirs(results_folder, exist_ok=True)
    
    ps_values = range(8, 33, 2)  # patch size from 10 to 20
    iterations = 100
    ite_per_transfer_values = [1, 2, 5, 10, 20, 25, 50, 100]  # ite_per_transfer from 1 to 100

    for d_value in d_values:
        results_file = os.path.join(results_folder, f"test_results_d_{d_value}.json")

        if os.path.exists(results_file):
            with open(results_file, 'r') as file:
                all_results = json.load(file)
            print(f"Loaded existing results from file for dimension {d_value}.")
        else:
            all_results = []
            for ps in ps_values:
                ps_results = []
                for ite_per_transfer in ite_per_transfer_values:
                    print(f"Running test for patch size {ps}, dimension {d_value}, ite_per_transfer {ite_per_transfer}...")
                    output = run_test(ps, d_value, ite_per_transfer)
                    result = parse_results(output)
                    ps_results.append(result[0] if result else {})
                all_results.append(ps_results)
            
            with open(results_file, 'w') as file:
                json.dump(all_results, file)
            print(f"Test results saved to file for dimension {d_value}.")

        # Plot the results for each ite_per_transfer value
        plot_results(all_results, ps_values, ite_per_transfer_values, d_value, results_folder)

if __name__ == "__main__":
    main()
