require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')

describe AlleleFrequency do
  before(:each) do
    @valid_attributes = {
      :freq => ,
      :population_number => 1
    }
  end

  it "should create a new instance given valid attributes" do
    AlleleFrequency.create!(@valid_attributes)
  end
end
