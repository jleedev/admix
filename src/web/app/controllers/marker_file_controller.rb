class MarkerFileController < ApplicationController
  def index
    @files = LocusFile.find :all, :order => 'updated_at DESC'
  end
end
