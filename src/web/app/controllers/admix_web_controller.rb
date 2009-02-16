class AdmixWebController < ApplicationController
  def index
    @title = "Admix web"
    @genotype_files = GenotypeFile.find :all, :order => 'updated_at DESC'
  end

  def create
  end
end
