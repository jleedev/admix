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
      catch :no_results do
        (yield :loc => loc, :ped => ped).should == out
      end
    end
  end

  describe Wrapper do
    it_should_behave_like "make test data available"

    it "should pass the admix test" do
      run_test do |data|
        Wrapper.call data
      end
    end

    it "should complain on empty input" do
      lambda {
        Wrapper.call :loc => "", :ped => ""
      }.should raise_error AdmixError
    end
  end

  describe AdmixWeb do
    it_should_behave_like "make test data available"

    it "should have a / page" do
      get_it "/"
      @response.should be_ok
    end

    it "should also pass the admix test" do
      run_test do |data|
        data.update :type => :text
        post_it "/admix", data
        @response.content_type.should == 'text/plain'
        @response.body
      end
    end

    it "should return some html" do
      run_test do |data|
        data.update :type => :html
        post_it "/admix", data
        @response.content_type.should == 'text/html'
        throw :no_results
      end
    end

    it "should report errors" do
      post_it "/admix"
      @response.should_not be_ok
    end

    it "'s prettify method should work" do
      (AdmixWeb::prettify "A   B  C\n1  2 3\n").should == "<tr><td>A<td>B<td>C\n<tr><td>1<td>2<td>3"
    end
  end

end
