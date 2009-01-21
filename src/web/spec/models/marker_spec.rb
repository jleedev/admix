require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')

describe Marker do
  before(:each) do
    @valid_attributes = {
      :name => 
    }
  end

  it "should create a new instance given valid attributes" do
    Marker.create!(@valid_attributes)
  end
end
