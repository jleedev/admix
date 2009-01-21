require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')

describe LocusFile do
  before(:each) do
    @valid_attributes = {
      :population => 1
    }
  end

  it "should create a new instance given valid attributes" do
    LocusFile.create!(@valid_attributes)
  end
end
