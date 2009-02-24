class Job < ActiveRecord::Base
  belongs_to :locus_file
  belongs_to :genotype_file
end
