$: << '../admixture/lib'
require 'admix'

class Job < ActiveRecord::Base
  belongs_to :locus_file
  belongs_to :genotype_file
  validates_presence_of :locus_file, :genotype_file

  def execute
    results = Admix::Wrapper.call :loc => locus_file.export, :ped => genotype_file.data
  end

  def name
    "#{locus_file.name} #{genotype_file.name}"
  end
end
