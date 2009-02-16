class GenotypeFile < ActiveRecord::Base
  belongs_to :user
  validates_presence_of :name,:data
end
