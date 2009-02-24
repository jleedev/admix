class Job < ActiveRecord::Base
  belongs_to :locus_file
  belongs_to :genotype_file
  validates_presence_of :locus_file, :genotype_file
end
