require 'rubygems'

require 'admix'

ENV['PATH'] = File.join(File.dirname(__FILE__), %w(.. bin)) + ':' + ENV['PATH']

def run_test
  loc,ped,out = %w(loc ped out).collect do |ext|
    File.read "orig/Test/admix-test." + ext
  end
  catch :no_results do
    (yield :loc => loc, :ped => ped).should == out
  end
end

describe Admix::Wrapper do
  it "should pass the admix test" do
    run_test do |data|
      Admix::Wrapper.call data
    end
  end

  it "should complain on empty input" do
    lambda {
      Admix::Wrapper.call :loc => "", :ped => ""
    }.should raise_error Admix::AdmixError
  end
end
