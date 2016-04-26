"""
This file contains helper functions used by many of the experiment scripts.
"""

from optparse import OptionParser
import os
from os import path
import subprocess
import sys

import plot_continuous_monitor

# Copy a file from a given host through scp, throwing an exception if scp fails.
def scp_from(host, identity_file, username, remote_file, local_file):
  subprocess.check_call("scp -q -o StrictHostKeyChecking=no -i {} '{}@{}:{}' '{}'".format(
    identity_file, username, host, remote_file, local_file), shell=True)

def ssh_get_stdout(host, identity_file, username, command):
  if "ec2" in host:
    command = "source /root/.bash_profile; {}".format(command)
  ssh_command = "ssh -t -o StrictHostKeyChecking=no -i {} {}@{} '{}'".format(
    identity_file, username, host, command)
  return subprocess.Popen(ssh_command, stdout=subprocess.PIPE, shell=True).communicate()[0]

def copy_latest_continuous_monitor(hostname, identity_file, filename_prefix, username):
  """ Copies logs back from a Spark cluster.

  This script copies the JSON event log and JSON continuous monitor back from a Spark
  driver and Spark executor, respectively, to the local machine.  Returns a two-item
  tuple with the name of the event log and the continuous monitor log.
  """
  continuous_monitor_relative_filename = ssh_get_stdout(
    hostname,
    identity_file,
    username,
    "ls -t /tmp/ | grep continuous_monitor | head -n 1").strip("\n").strip("\r")
  continuous_monitor_filename = "/tmp/{}".format(continuous_monitor_relative_filename)
  local_continuous_monitor_file = "{}_executor_monitor".format(filename_prefix)
  print "Copying continuous monitor from file {} on host {} back to {}".format(
    continuous_monitor_filename, hostname, local_continuous_monitor_file)
  scp_from(
    hostname,
    identity_file,
    username,
    continuous_monitor_filename,
    local_continuous_monitor_file)

  print "Plotting continuous monitor"
  plot_continuous_monitor.plot_continuous_monitor(local_continuous_monitor_file, open_graphs=True)
  return local_continuous_monitor_file

def plot_continuous_monitors(log_dir):
  """ Plots all of the continuous monitors in the provided directory. """
  for log_filename in os.listdir(log_dir):
    if log_filename.endswith("executor_monitor"):
      plot_continuous_monitor.plot_continuous_monitor(
        path.join(log_dir, log_filename), use_gnuplot=True)
