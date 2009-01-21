class Marker < ActiveRecord::Base
  belongs_to :locus_file
  has_many :allele_frequency
end
