import argparse
import inspect
import json
import os
import subprocess
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

BYTES_PER_GIGABYTE = float(1024 * 1024 * 1024)

def continuous_monitor_col(key, continuous_monitor):
  """
  For a given key, returns a list of data from continuous monitor entries
  corresponding to that key.
  """
  return [data[key] for data in continuous_monitor]

def plot_continuous_monitor(filename, open_graphs=False, use_gnuplot=False):
  continuous_monitor_data = []
  out_filename = "%s_utilization" % filename
  if use_gnuplot:
    out_file = open(out_filename, "w")

  # Get the location of the monotasks-scripts repository by getting the directory containing the
  # file that is currently being executed.
  scripts_dir = os.path.dirname(inspect.stack()[0][1])

  start = -1
  at_beginning = True
  for (i, line) in enumerate(open(filename, "r")):
    try:
      json_data = json.loads(line)
    except ValueError:
      # This typically happens at the end of the file, which can get cutoff when the job stops.
      print "Stopping parsing due to incomplete line"
      if not at_beginning:
        break
      else:
        # There are some non-JSON lines at the beginning of the file.
        print "Skipping non-JSON line at beginning of file: %s" % line
        continue
    at_beginning = False
    time = json_data["Current Time"]
    if start == -1:
      start = time
    disk_utilizations = json_data["Disk Utilization"]["Device Name To Utilization"]
    xvdf_utilization = get_util_for_disk(disk_utilizations, "xvdf")
    xvdb_utilization= get_util_for_disk(disk_utilizations, "xvdb")
    xvdf_total_utilization = xvdf_utilization["Disk Utilization"]
    xvdb_total_utilization = xvdb_utilization["Disk Utilization"]
    xvdf_read_throughput = xvdf_utilization["Read Throughput"]
    xvdf_write_throughput = xvdf_utilization["Write Throughput"]
    xvdb_read_throughput = xvdb_utilization["Read Throughput"]
    xvdb_write_throughput = xvdb_utilization["Write Throughput"]
    cpu_utilization = json_data["Cpu Utilization"]
    cpu_system = cpu_utilization["Total System Utilization"]
    cpu_total = (cpu_utilization["Total User Utilization"] +
      cpu_utilization["Total System Utilization"])
    network_utilization = json_data["Network Utilization"]
    bytes_received = network_utilization["Bytes Received Per Second"]
    running_compute_monotasks = 0
    if "Running Compute Monotasks" in json_data:
      running_compute_monotasks = json_data["Running Compute Monotasks"]

    xvdf_running_disk_monotasks = 0
    xvdb_running_disk_monotasks = 0
    if "Running Disk Monotasks" in json_data:
      # Parse the number of currently running disk monotasks for each disk.
      for running_disk_monotasks_info in json_data["Running Disk Monotasks"]:
        running_disk_monotasks = running_disk_monotasks_info["Running And Queued Monotasks"]
        disk_name = running_disk_monotasks_info["Disk Name"]
        if "xvdf" in disk_name:
          xvdf_running_disk_monotasks = running_disk_monotasks
        elif "xvdb" in disk_name:
          xvdb_running_disk_monotasks = running_disk_monotasks

    running_macrotasks = 0
    if "Running Macrotasks" in json_data:
      running_macrotasks = json_data["Running Macrotasks"]
    local_running_macrotasks = 0
    if "Local Running Macrotasks" in json_data:
      local_running_macrotasks = json_data["Local Running Macrotasks"]
    gc_fraction = 0
    if "Fraction GC Time" in json_data:
      gc_fraction = json_data["Fraction GC Time"]
    outstanding_network_bytes = 0
    if "Outstanding Network Bytes" in json_data:
      outstanding_network_bytes = json_data["Outstanding Network Bytes"]
    if bytes_received == "NaN" or bytes_received == "Infinity":
      continue
    bytes_transmitted = network_utilization["Bytes Transmitted Per Second"]
    if bytes_transmitted == "NaN" or bytes_transmitted == "Infinity":
      continue
    if str(cpu_total).find("NaN") > -1 or str(cpu_total).find("Infinity") > -1:
      continue
    macrotasks_in_network = 0
    if "Macrotasks In Network" in json_data:
      macrotasks_in_network = json_data["Macrotasks In Network"]
    macrotasks_in_compute = 0
    if "Macrotasks In Compute" in json_data:
      macrotasks_in_compute = json_data["Macrotasks In Compute"]
    macrotasks_in_disk = 0
    if "Macrotasks In Disk" in json_data:
      macrotasks_in_disk = json_data["Macrotasks In Disk"]
    free_heap_memory = 0
    if "Free Heap Memory Bytes" in json_data:
      free_heap_memory = json_data["Free Heap Memory Bytes"]
    free_off_heap_memory = 0
    if "Free Off-Heap Memory Bytes" in json_data:
      free_off_heap_memory = json_data["Free Off-Heap Memory Bytes"]

    data = {
      'time': time - start,
      'xvdf utilization': xvdf_total_utilization,
      'xvdb utilization': xvdb_total_utilization,
      'cpu utilization': cpu_total / 8.0, # Cores in machine?
      'bytes received': bytes_received / 125000000.,
      'bytes transmitted': bytes_transmitted / 125000000., # what is the 125000000 magic number?
      'running compute monotasks': running_compute_monotasks,
      'running monotasks': running_macrotasks,
      'gc fraction': gc_fraction,
      'outstanding network bytes': outstanding_network_bytes / (1024 * 1024),
      'macrotasks in network': macrotasks_in_network,
      'macrotasks in compute': macrotasks_in_compute,
      'cpu system': cpu_system / 8.0,
      'macrotasks in disk': macrotasks_in_disk,
      'xvdf read throughput': xvdf_read_throughput,
      'xvdf write throughput': xvdf_write_throughput,
      'xvdb read throughput': xvdb_read_throughput,
      'xvdb write throughput': xvdb_write_throughput,
      'xvdf running disk monotasks': xvdf_running_disk_monotasks,
      'xvdb running disk monotasks': xvdb_running_disk_monotasks,
      'free heap memory': free_heap_memory / BYTES_PER_GIGABYTE,
      'free off heap memory': free_off_heap_memory / BYTES_PER_GIGABYTE
    }
    continuous_monitor_data.append(data)
    if use_gnuplot:
      write_data(out_file, gnuplot_data_line(data))

  if use_gnuplot:
    out_file.close()
    plot_gnuplot(filename, open_graphs)
  else:
    plot_matplotlib(continuous_monitor_data, filename, open_graphs)

def get_util_for_disk(disk_utils, disk):
  """
  Returns the disk utilization metrics for the specified disk, given the utilization information
  for all disks, or None if the desired disk cannot be found.
  """
  for disk_util in disk_utils:
    if disk in disk_util:
      return disk_util[disk]
  return None

def gnuplot_data_line(data_hash):
  """ Helper function to order continuous montior data for gnuplot to use """
  return [data_hash['time'],
          data_hash['xvdf utilization'],
          data_hash['xvdb utilization'],
          data_hash['cpu utilization'],
          data_hash['bytes received'],
          data_hash['bytes transmitted'],
          data_hash['running compute monotasks'],
          data_hash['running monotasks'],
          data_hash['gc fraction'],
          data_hash['outstanding network bytes'],
          data_hash['macrotasks in network'],
          data_hash['macrotasks in compute'],
          data_hash['cpu system'],
          data_hash['macrotasks in disk'],
          data_hash['xvdf read throughput'],
          data_hash['xvdf write throughput'],
          data_hash['xvdb read throughput'],
          data_hash['xvdb write throughput'],
          data_hash['xvdf running disk monotasks'],
          data_hash['xvdb running disk monotasks'],
          data_hash['free heap memory'],
          data_hash['free off heap memory']]

def plot_gnuplot(file_prefix, open_graphs):
  """
  Creates modified copies of gnuplot files, and associated pdf graphs for several
  attributes from the continuous monitor.
  """
  # Get the location of the monotasks-scripts repository by getting the directory containing the
  # file that is currently being executed.
  scripts_dir = os.path.dirname(inspect.stack()[0][1])
  attributes = ['utilization', 'disk_utilization', 'monotasks', 'memory']

  for attribute in attributes:
    plot_gnuplot_attribute(file_prefix, attribute, scripts_dir)

  if open_graphs:
    for attribute in attributes:
      subprocess.check_call('open {0}_{1}.pdf'.format(file_prefix, attribute), shell=True)

def plot_gnuplot_attribute(file_prefix, attribute, scripts_dir):
  """ Create a gnuplot file and associated pdf for a some attribute (like disk_utilization) """
  data_filename = '{0}_utilization'.format(file_prefix)
  plot_filename = '{0}_{1}.gp'.format(file_prefix, attribute)
  pdf_filename = '{0}_{1}.pdf'.format(file_prefix, attribute)
  plot_file = open(plot_filename, 'w')

  for line in open(os.path.join(scripts_dir, 'gnuplot_files/plot_{0}_base.gp'.format(attribute)), 'r'):
    new_line = line.replace('__OUT_FILENAME__', pdf_filename).replace('__NAME__', data_filename)
    plot_file.write(new_line)
  plot_file.close()

  subprocess.check_call('gnuplot {0}'.format(plot_filename), shell=True)

def plot_matplotlib(cm_data, file_prefix, open_graphs):
  disk_utilization_params = ['xvdf utilization',
                              'xvdb utilization',
                              'xvdf write throughput',
                              'xvdf read throughput',
                              'xvdb write throughput',
                              'xvdb read throughput']
  memory_params = ['free heap memory',
                    'free off heap memory']
  monotasks_params = ['macrotasks in network',
                      'macrotasks in compute',
                      'macrotasks in disk',
                      'running monotasks',
                      'gc fraction',
                      'outstanding network bytes']
  utilization_params = ['cpu utilization',
                        'xvdf utilization',
                        'xvdb utilization',
                        'cpu system',
                        'gc fraction']
  xvdf_params = ['xvdf running disk monotasks',
                 'xvdf write throughput',
                 'xvdf read throughput',
                 'xvdf utilization']
  xvdb_params = ['xvdb running disk monotasks',
                 'xvdb write throughput',
                 'xvdb read throughput',
                 'xvdb utilization']

  def plot_params(params_to_plot, title):
    """
    Creates a matplotlib graph using continuous monitor data.Time is the x axis and data corresponding to
    each parameter is used to generate a new line on the line graph.
    """
    filename = '{0}_{1}.pdf'.format(file_prefix, '_'.join(title.lower().split()))
    handles = []
    time = continuous_monitor_col('time', cm_data)
    for key in params_to_plot:
      handle, = plt.plot(time, continuous_monitor_col(key, cm_data), label=key)
      handles.append(handle)
    plt.legend(handles = handles)
    plt.title(title)
    pdf.savefig()
    if open_graphs:
      plt.show()
    else:
      plt.close()

  with PdfPages('{0}_graphs.pdf'.format(file_prefix)) as pdf:
    plot_params(disk_utilization_params,
                'Disk Utilization')
    plot_params(memory_params,
                'Memory')
    plot_params(monotasks_params,
                'Monotasks')
    plot_params(utilization_params,
                'Utilization')
    plot_params(xvdb_params,
                'xvdd Utilization')
    plot_params(xvdf_params,
                'xvdf Utilization')

def plot_single_disk(filename, utilization_filename, disk_to_plot, disks_to_skip, scripts_dir):
  """
  Plots the utilization for a single disk, ignoring the utilization for any disks in
  disks_to_skip.
  """
  disk_plot_filename_prefix = "%s_%s_disk_utilization" % (filename, disk_to_plot)
  disk_plot_filename = "%s.gp" % disk_plot_filename_prefix
  disk_plot_output = "%s.pdf" % disk_plot_filename_prefix
  disk_plot_file = open(disk_plot_filename, "w")
  for line in open(os.path.join(scripts_dir, "gnuplot_files/plot_disk_utilization_base.gp"), "r"):
    skip = False
    for disk_to_skip in disks_to_skip:
      if line.find(disk_to_skip) != -1:
        skip = True
    if not skip:
      new_line = line.replace("__OUT_FILENAME__", disk_plot_output).replace(
        "__NAME__", utilization_filename)
      disk_plot_file.write(new_line)

  disk_plot_file.close()
  subprocess.check_call("gnuplot %s" % disk_plot_filename, shell=True)
  return disk_plot_output

def write_data(out_file, data):
  stringified = [str(x) for x in data]
  out_file.write("\t".join(stringified))
  out_file.write("\n")

def parse_args():
  parser = argparse.ArgumentParser(description="Plots Spark continuous monitor logs.")
  parser.add_argument(
    "-f", "--filename", help="The path to a continuous monitor log file.", required=True)
  parser.add_argument(
    "-o", "--open_graphs", help="open generated graphs", action="store_true")
  parser.add_argument(
    "-g", "--gnuplot", help="generate graphs with gnuplot", action="store_true")
  parser.set_defaults(gnuplot=False, open_graphs=False)
  return parser.parse_args()

def main():
  args = parse_args()
  plot_continuous_monitor(args.filename, args.open_graphs, args.gnuplot)

if __name__ == "__main__":
  main()
