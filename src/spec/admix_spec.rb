require 'admix'

describe Admix::AdmixWrapper do
  it "should pass the admix tests" do
    loc,ped,out = %w(loc ped out).each do |ext|
      File.read "../orig/Test/admix-test." + ext
    end
    @admix = Admix::AdmixWrapper.new :loc => loc, :ped => ped
    @admix.out.should == out
  end
end
