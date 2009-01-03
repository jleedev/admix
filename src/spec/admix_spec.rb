require 'admix'

module Admix

  describe AdmixWrapper do
    it "should pass the admix tests" do
      loc,ped,out = %w(loc ped out).collect do |ext|
        File.read "orig/Test/admix-test." + ext
      end
      admix = AdmixWrapper.new :loc => loc, :ped => ped
      admix.out.should == out
    end

    it "should complain on empty input" do
      lambda {
        AdmixWrapper.new :loc => "", :pred => ""
      }.should raise_error AdmixError
    end
  end

end
