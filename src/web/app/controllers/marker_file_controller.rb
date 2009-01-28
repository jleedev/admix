class MarkerFileController < ApplicationController
  def index
    @title = "Marker files"
    @files = LocusFile.find :all, :order => 'updated_at DESC'
  end

  def create
    @title = "Create marker file"
  end
end
