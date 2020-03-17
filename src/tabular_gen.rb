#!/usr/bin/env ruby

PREFIX = "logs/"
MIXING_PREFIX = "../../mixing/logs/"
SUFFIX = ["_0.1k", "_0.2k", "_0.3k", "_full"]
GSET_RNG = 1..67
REGEX = /iter: (\d+), cut: (\d+\.\d+e[+-]\d+), time: (\d+\.\d+)/
REGEX_MIXING = /iter: (\d+), cut: (\d+\.\d+), time: (\d+\.\d+)/
K_CAP = 200

def get_csr_size(fname)
  File.open(fname).each_with_index do |line, i|
    if (i == 2)
      return line.split(/\s+/)[1].to_i
    end
  end
end

GSET_RNG.each do |i|
  n = get_csr_size("../data/G/G#{i}")
  kk = [(Math::sqrt(2 * n)).ceil, K_CAP].min
  str = sprintf("G%d & %5d & %5d ", i, n, kk)
  # rbr
  SUFFIX.each do |s|
    fname = PREFIX + "G#{i}#{s}.log"
    full_text = File.open(fname, "r").read
    REGEX =~ full_text
    str += sprintf("& %6d & %13.6e & %6.2f ", $1.to_i, $2.to_f, $3.to_f)
  end
  # mixing method
  fname = MIXING_PREFIX + "G#{i}.log"
  full_text = File.open(fname, "r").read
  REGEX_MIXING =~ full_text
  str += sprintf("& %6d & %13.6e & %6.2f ", $1.to_i, $2.to_f, $3.to_f)
  str += "\\\\ \n"
  print str
end
