class AdmixWebController < ApplicationController
  def index
    @title = "Admix web"
    @genotype_files = GenotypeFile.find :all, :order => 'updated_at DESC'
  end

  def create
    @title = "Upload genotype file"
    @genotype_file = GenotypeFile.new params[:genotype_file]
    @genotype_file.data = @genotype_file.data.read if @genotype_file.data
    if request.post? and @genotype_file.save
      redirect_to :action => :index
    end
  end

  def show
    @genotype_file = GenotypeFile.find_by_id params[:id]
    @title = "Genotype file \"#{@genotype_file.name}\""
  end

  def create_job
    @title = "Create job"
    @job = Job.new params[:job]
    if request.post? and @job.save
      redirect_to :action => :index
    end
  end
end
