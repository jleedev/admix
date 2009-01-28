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

  def show
    @locus_file = LocusFile.find_by_id params[:id]
    @title = "Marker file \"#{@locus_file.name}\""
  end

  def create_marker
    @locus_file = LocusFile.find_by_id params[:id]
    marker = Marker.new params[:marker]
    marker.locus_file = @locus_file
    marker.save
    params[:allele_frequency].each do |k,v|
      f = AlleleFrequency.new :marker => marker, :population_number => k.to_i, :freq => v[:freq]
      f.save
    end
    redirect_to :action => :show, :id => @locus_file.id
  end
end
