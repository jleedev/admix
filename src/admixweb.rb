$: << File.join(File.expand_path(File.dirname(__FILE__)), '..', 'lib')

require 'rubygems'
require 'sinatra'

require 'admix'

module Admix

  get '/' do <<EOHTML
Let's mix up some ads.
EOHTML
  end

  post '/admix' do |*args|
    admix = AdmixWrapper.new :loc => params[:loc], :ped => params[:ped]
    admix.out
  end

end
