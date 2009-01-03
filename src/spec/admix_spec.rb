$: << File.join(File.expand_path(File.dirname(__FILE__)), '..', 'lib')

require 'rubygems'
require 'sinatra'
require 'sinatra/test/rspec'

require 'admix'
require 'admixweb'

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

  describe "AdmixWeb" do
    it "should have a / page" do
      get_it "/"
      @response.should be_ok
    end

    it "should also pass the admix tests" do
      loc,ped,out = %w(loc ped out).collect do |ext|
        File.read "orig/Test/admix-test." + ext
      end
      post_it "/admix", :loc => loc, :ped => ped
      @response.body.should == out
    end
  end

end
