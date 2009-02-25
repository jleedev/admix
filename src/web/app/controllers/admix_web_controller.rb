$: << '../admixture/lib'
require 'admix'

class AdmixWebController < ApplicationController
  def index
    @title = "Admix web"
    @genotype_files = GenotypeFile.find :all, :order => 'updated_at DESC'
    @jobs = Job.find :all, :order => 'updated_at DESC'
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
    if request.post?
      locus_file = LocusFile.find_by_id params[:job][:locus_file_id] if params[:job]
      genotype_file = GenotypeFile.find_by_id params[:job][:genotype_file_id] if params[:job]
      @job = Job.new params[:job]
      @job.locus_file_id = locus_file
      @job.genotype_file_id = genotype_file
      begin
        @job.execute
      rescue Admix::AdmixError => e
        @job.results = e.message
      end
      if @job.save
        redirect_to :action => :index
      end
    end
  end

  def show_job
    @job = Job.find_by_id params[:id]
    @title = "Job \"#{@job.name}\""
  end
end
