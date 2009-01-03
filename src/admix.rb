require 'tempfile'

module Admix
  class AdmixWrapper
    def initialize args
      Tempfile.open 'admix.loc' do |loc|
        (loc << args[:loc]).flush
        Tempfile.open 'admix.ped' do |ped|
          (ped << args[:ped]).flush
          Tempfile.open 'admix.out' do |out|
            result = `orig/admix -q -g .1 -M "0" #{loc.path} #{ped.path} #{out.path} 2>&1`
            case $?.exitstatus
            when 0
              return out.read
            when 1
              throw result
            end
          end
        end
      end

      IO.popen "orig/admix"
    end

    def out
      return @out
    end
  end
end
