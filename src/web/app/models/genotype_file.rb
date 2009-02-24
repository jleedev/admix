class GenotypeFile < ActiveRecord::Base
  belongs_to :user
  has_many :job
  validates_presence_of :name,:data
end
