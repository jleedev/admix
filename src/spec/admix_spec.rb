require 'rubygems'
require 'sinatra'
require 'sinatra/test/rspec'

require 'admix'
require 'admixweb'

module Admix

  describe "make test data available", :shared => true do
    def run_test
      loc,ped,out = %w(loc ped out).collect do |ext|
        File.read "orig/Test/admix-test." + ext
      end
      out.should == (yield :loc => loc, :ped => ped)
    end
  end

  describe AdmixWrapper do
    it_should_behave_like "make test data available"

    it "should pass the admix test" do
      run_test do |data|
        admix = AdmixWrapper.new data
        admix.out
      end
    end

    it "should complain on empty input" do
      lambda {
        AdmixWrapper.new :loc => "", :ped => ""
      }.should raise_error AdmixError
    end
  end

  describe "AdmixWeb" do
    it_should_behave_like "make test data available"

    it "should have a / page" do
      get_it "/"
      @response.should be_ok
    end

    it "should also pass the admix test" do
      run_test do |data|
        post_it "/admix", data
        @response.body
      end
    end

    it "should report errors" do
      post_it "/admix"
      @response.should_not be_ok
    end
  end

end
