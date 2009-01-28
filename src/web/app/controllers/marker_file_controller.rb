class MarkerFileController < ApplicationController
  def index
    @title = "Marker files"
    @files = LocusFile.find :all, :order => 'updated_at DESC'
  end

  def create
    @title = "Create marker file"
    @locus_file = LocusFile.new params[:locus_file]
    if request.post? and @locus_file.save
      redirect_to :action => :index
    end
  end
end
