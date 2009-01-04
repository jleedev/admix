require 'tempfile'

module Admix

  # A simple wrapper around the command line admix utility. This takes
  # strings, passes them to the admix command, and then returns its
  # output as a string.
  class Wrapper
    def self.call args
      Tempfile.open 'admix.loc' do |loc|
      (loc << args[:loc]).flush
      Tempfile.open 'admix.ped' do |ped|
      (ped << args[:ped]).flush
      Tempfile.open 'admix.out' do |out|
        result = `orig/admix -q -g .1 -M "0" #{loc.path} #{ped.path} #{out.path} 2>&1`
        case $?.exitstatus
        when 0
          return out.read
        else
          raise AdmixError::new result.strip
        end
      end
      end
      end
    end
  end

  class AdmixError < Exception; end

end
